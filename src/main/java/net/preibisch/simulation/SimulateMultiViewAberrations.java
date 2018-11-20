/*-
 * #%L
 * Library for simulating a multi-view acquisition including
 * attenuation, convolution, reduced sampling and poission noise.
 * %%
 * Copyright (C) 2014 - 2017 Multiview Simulation developers.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
package net.preibisch.simulation;

import java.util.ArrayList;
import java.util.Date;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import ij.IJ;
import ij.ImageJ;
import mpicbg.models.AffineModel3D;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.gradient.HessianMatrix;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.multithreading.SimpleMultiThreading;
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory;
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory.Boundary;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.Util;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

/**
 * Code for simulating a multi-view acquisition including attenuation, convolution, reduced sampling and poisson noise
 * 
 * This software is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see http://www.gnu.org/licenses/.
 * 
 * @author Stephan Preibisch (stephan.preibisch@gmx.de)
 */
public class SimulateMultiViewAberrations
{
	final static Random rnd = new Random( 464232194 );
	final public static float minValue = 0.0001f;
	public final static float avgIntensity = 1;

	public static AffineModel3D axisRotation( final Interval in, final int axis, final int degrees )
	{
		// translate so that the center of the image is 0,0,0
		final AffineModel3D translate1 = new AffineModel3D();
		translate1.set( 1, 0, 0, -( in.max( 0 ) - in.min( 0 ) )/2,
					    0, 1, 0, -( in.max( 1 ) - in.min( 1 ) )/2,
					    0, 0, 1, -( in.max( 2 ) - in.min( 2 ) )/2 );
		
		// rotate around an axis
		final AffineModel3D rot = new AffineModel3D();
		rot.rotate( axis, (float)Math.toRadians( degrees ) );

		// translate back to the center
		final AffineModel3D translate2 = new AffineModel3D();
		translate2.set( 1, 0, 0, ( in.max( 0 ) - in.min( 0 ) )/2,
					    0, 1, 0, ( in.max( 1 ) - in.min( 1 ) )/2,
					    0, 0, 1, ( in.max( 2 ) - in.min( 2 ) )/2 );

		translate1.preConcatenate( rot );
		translate1.preConcatenate( translate2 );

		return translate1;
	}

	public static Img< FloatType > rotateAroundAxis( final RandomAccessibleInterval< FloatType > in, final int axis, final int degrees )
	{
		// final model
		final AffineModel3D affine = axisRotation( in, axis, degrees ).createInverse();

		final Img< FloatType > out = new ArrayImgFactory< FloatType >().create( in, new FloatType() );

		final NLinearInterpolatorFactory< FloatType > factory = new NLinearInterpolatorFactory< FloatType >();
		final RealRandomAccessible< FloatType > interpolant = Views.interpolate( Views.extendZero( in ), factory );
		final RealRandomAccess< FloatType > realRandomAccess = interpolant.realRandomAccess();

		final Cursor< FloatType > c = out.localizingCursor();
		final int[] l = new int[ out.numDimensions() ];
		final double[] lf = new double[ out.numDimensions() ];
		
		while ( c.hasNext() )
		{
			c.fwd();
			c.localize( l );

			lf[ 0 ] = l[ 0 ];
			lf[ 1 ] = l[ 1 ];
			lf[ 2 ] = l[ 2 ];
			
			affine.applyInPlace( lf );
			realRandomAccess.setPosition( lf );
			
			c.get().set( realRandomAccess.get() );
		}

		return out;
	}
	
	/**
	 * Scales the reduced lightsheet acquisition back to isotropic size
	 * 
	 * @param randomAccessible - the input
	 * @param inc - every n'th 
	 * @return - the isotropic image
	 */
	public static Img< FloatType > makeIsotropic( final RandomAccessibleInterval< FloatType > randomAccessible, final int inc )
	{
		final long[] dim = new long[]{ randomAccessible.dimension( 0 ), randomAccessible.dimension( 1 ), ( randomAccessible.dimension( 2 ) - 1 ) * inc + 1 };
		final Img< FloatType > img = new ArrayImgFactory< FloatType >().create( dim, new FloatType() );
		
		final NLinearInterpolatorFactory< FloatType > factory = new NLinearInterpolatorFactory< FloatType >();
		final RealRandomAccessible< FloatType > interpolant = Views.interpolate( Views.extendMirrorSingle( randomAccessible ), factory );
		final RealRandomAccess< FloatType > realRandomAccess = interpolant.realRandomAccess();
		
		final Cursor< FloatType > c = img.localizingCursor();
		final int[] l = new int[ img.numDimensions() ];
		final double[] lf = new double[ img.numDimensions() ];
		
		while ( c.hasNext() )
		{
			c.fwd();
			c.localize( l );
			
			lf[ 0 ] = l[ 0 ];
			lf[ 1 ] = l[ 1 ];
			lf[ 2 ] = (float)l[ 2 ] / (float)inc;
			
			realRandomAccess.setPosition( lf );
			c.get().set( realRandomAccess.get() );
		}
		
		return img;
	}
	
	/**
	 * Scans the sample with a simulated lightsheet ... the width of the lightsheet is implicitly defined by the effective PSF
	 * 
	 * @param randomAccessible - input
	 * @param inc - every n'th 
	 * @param poissonSNR - which poisson SNR is desired?
	 * @return every n'th slice
	 */
	public static Img< FloatType > extractSlices( final RandomAccessibleInterval< FloatType > randomAccessible, final int inc, final float poissonSNR  )
	{
		return extractSlices( randomAccessible, inc, poissonSNR, rnd );
	}

	/**
	 * Scans the sample with a simulated lightsheet ... the width of the lightsheet is implicitly defined by the effective PSF
	 * 
	 * @param randomAccessible - input
	 * @param inc - every n'th 
	 * @param poissonSNR - which poisson SNR is desired?
	 * @param rnd - a random number generator
	 * @return every n'th slice
	 */
	public static Img< FloatType > extractSlices( final RandomAccessibleInterval< FloatType > randomAccessible, final int inc, final float poissonSNR, final Random rnd  )
	{
		final long[] dim = new long[]{ randomAccessible.dimension( 0 ), randomAccessible.dimension( 1 ), ( randomAccessible.dimension( 2 ) - 1 )/inc + 1 };
		final Img< FloatType > img = new ArrayImgFactory< FloatType >().create( dim, new FloatType() );

		final RandomAccess< FloatType > r = img.randomAccess();
		final int[] tmp = new int[ 3 ];

		IJ.showProgress( 0 );

		int countZ = 0;
		for ( int z = 0; z < randomAccessible.dimension( 2 ); z += inc )
		{
			final RandomAccessibleInterval< FloatType > slice = Views.hyperSlice( randomAccessible, 2, z );
			final Cursor< FloatType > c;
			
			if ( poissonSNR >= 0.0 )
				c = poissonProcess( slice, poissonSNR, rnd ).localizingCursor();
			else
				c = Views.iterable( slice ).localizingCursor();
			
			tmp[ 2 ] = countZ++;
			while ( c.hasNext() )
			{
				c.fwd();
				tmp[ 0 ] = c.getIntPosition( 0 );
				tmp[ 1 ] = c.getIntPosition( 1 );
				
				r.setPosition( tmp );
				r.get().set( c.get() );
			}

			IJ.showProgress( z, (int)randomAccessible.dimension( 2 ) );
		}
		
		return img;
	}
	
	public static Img< FloatType > poissonProcess( final RandomAccessibleInterval< FloatType > in, final float poissonSNR, final Random rnd )
	{
		final Img< FloatType > out = new ArrayImgFactory< FloatType >().create( in, new FloatType() );

		final Cursor< FloatType > c = out.localizingCursor();
		final RandomAccess< FloatType > r = in.randomAccess();
		
		while ( c.hasNext() )
		{
			c.fwd();
			r.setPosition( c );
			c.get().set( r.get() );
		}
		
		// based on an average intensity of 5 inside the sample
		Tools.poissonProcess( out, poissonSNR, rnd );
		
		return out;
	}
	
	public static Img< FloatType > convolve( final Img< FloatType > img, final Img< FloatType > psf, final ExecutorService service )
	{
		Tools.normImage( psf );
		final Img< FloatType > result = img.factory().create( img, img.firstElement() );
		final FFTConvolution< FloatType > conv = new FFTConvolution<FloatType>( img, psf, result, getFFTFactory( img ), service );
		
		// this fixes the wrong default kernel flipping in older versions of FFTConvolution 
		conv.setComputeComplexConjugate(false);
		conv.convolve();
		
		return result;
	}
	
	protected static ImgFactory< ComplexFloatType > getFFTFactory( final Img< ? extends RealType< ? > > img )
	{
		try
		{
			return img.factory().imgFactory( new ComplexFloatType() );
		}
		catch ( final IncompatibleTypeException e )
		{
			if ( img.size() > Integer.MAX_VALUE / 2 )
				return new CellImgFactory< ComplexFloatType >( 1024 );
			return new ArrayImgFactory< ComplexFloatType >();
		}
	}

	public static Img< FloatType > computeWeightImage( final RandomAccessibleInterval< FloatType > randomAccessible, final double delta )
	{
		// over which the cosine function spans
		final int cosineSpan = 40;
		
		// the weight image
		final Img< FloatType > img = new ArrayImgFactory< FloatType >().create( randomAccessible, new FloatType() );
		
		final Cursor< FloatType > c = img.localizingCursor();
		final int sizeY = (int)randomAccessible.dimension( 1 );
		
		while ( c.hasNext() )
		{
			c.fwd();
			
			int l = ( sizeY - c.getIntPosition( 1 ) - 1 );
			float value;
			
			if ( l < sizeY/2 )
			{
				value = 1.0f;
			}
			else if ( l > sizeY/2 + cosineSpan )
			{
				value = 0.0f;
			}
			else
			{
				final double pos = ( (double)(l - sizeY/2) / (double)cosineSpan ) * Math.PI;
				value = (float)( ( Math.cos( pos ) + 1.0 ) / 2.0 );
			}
			
			c.get().set( value );
		}
		
		return img;
	}
	
	public static Img< FloatType > attenuate3d( final RandomAccessibleInterval< FloatType > randomAccessible, final double delta )
	{
		// the attenuated image
		final Img< FloatType > img = new ArrayImgFactory< FloatType >().create( randomAccessible, new FloatType() );
		
		// make a plane that only contains x & z, this is the origin for the attenuation along y
		final RandomAccessibleInterval< FloatType > startMatrix = Views.hyperSlice( randomAccessible, 1, 0 );
		final Cursor< FloatType > c = Views.iterable( startMatrix ).localizingCursor();
		
		final RandomAccess< FloatType > rIn = randomAccessible.randomAccess();
		final RandomAccess< FloatType > rOut = img.randomAccess();
		final int[] l = new int[ 3 ];
		
		while ( c.hasNext() )
		{
			c.fwd();
			
			l[ 0 ] = c.getIntPosition( 0 );
			l[ 1 ] = (int)randomAccessible.dimension( 1 ) - 1;
			l[ 2 ] = c.getIntPosition( 1 );
			
			rIn.setPosition( l );
			rOut.setPosition( l );
			
			// light intensity goes down from 1...0 depending on the image intensities
			double n = 1;
			
			for ( int y = 0; y < randomAccessible.dimension( 0 ); ++y )
			{
				final double v = rIn.get().get();
				
				// probability that light is absorbed at this point
				final double phiN = v * delta * n;
				
				// n is >=0
				n = Math.max( n - phiN, 0 );
				
				// set attenuated intensity
				rOut.get().set( (float) ( v * n ) );
				
				rIn.bck( 1 );
				rOut.bck( 1 );
			}
		}
		
		return img;
	}

	public static void drawSimpleImage( final RandomAccessibleInterval< FloatType > randomAccessible )
	{
		final Cursor< FloatType > c = Views.iterable( randomAccessible ).localizingCursor();

		while( c.hasNext() )
		{
			c.fwd();
			boolean all = true;
			//for ( int d = 0; d < n; ++d )
			//	if ( c.getIntPosition( d ) > 100 || c.getIntPosition( d ) < 10 )
			//		all = false;

			if ( c.getIntPosition( 0 ) < randomAccessible.dimension( 0 ) / 2 )
			{
				if ( c.getIntPosition( 0 ) > c.getIntPosition( 1 ) )
					all = true;
				else
					all = false;
			}
			else
			{
				if ( randomAccessible.dimension( 0 ) - c.getIntPosition( 0 ) > c.getIntPosition( 1 ) )
					all = true;
				else
					all = false;
			}

			if ( all )
				c.get().set( 1 );
			else
				c.get().set( 0 );
		}
	}

	public static Img< FloatType > refract3d( final RandomAccessibleInterval< FloatType > randomAccessible )
	{
		// the refracted image
		final Img< FloatType > img = new ArrayImgFactory< FloatType >( new FloatType() ).create( randomAccessible );
		
		// make a plane that only contains x & z, this is the origin for the attenuation along y
		final RandomAccessibleInterval< FloatType > startMatrix = Views.hyperSlice( randomAccessible, 1, 0 );
		final Cursor< FloatType > c = Views.iterable( startMatrix ).localizingCursor();
		
		final RandomAccess< FloatType > rIn = Views.extendMirrorSingle( randomAccessible ).randomAccess();
		final RandomAccess< FloatType > rOut = img.randomAccess();
		final int[] l = new int[ 3 ];

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		while ( c.hasNext() )
		{
			c.fwd();

			l[ 0 ] = c.getIntPosition( 0 );
			l[ 1 ] = (int)randomAccessible.dimension( 1 ) - 1;
			l[ 2 ] = c.getIntPosition( 1 );

			if ( l[ 0 ] % 13 != 0 )
				continue;

			rIn.setPosition( l );
			rOut.setPosition( l );

			double vNext;

			for ( int y = 1; y < randomAccessible.dimension( 0 ); ++y )
			{
				// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)

				Hessian.computeHessianMatrix3D( rIn, matrix );

				double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

				double nx = eigenVector[ 0 ];
				double ny = eigenVector[ 1 ];
				double nz = eigenVector[ 2 ];

				// current direction of the ray
				final double bx = 0;
				final double by = 1;
				final double bz = 0;

				double dotP = Math.acos( ( nx*bx + ny*by + nz*bz) / ( Math.sqrt( nx*nx + ny*ny + nz*nz ) * Math.sqrt( bx*bx + by*by + bz*bz ) ) );

				if ( dotP >= Math.PI / 2 )
				{
					// invert normal vector & eigenvalue
					nx *= -1;
					ny *= -1;
					nz *= -1;
					ev *= -1;

					// adjust angle
					dotP -= Math.PI / 2;
					//dotP = Math.acos( ( nx*bx + ny*by + nz*bz) / ( Math.sqrt( nx*nx + ny*ny + nz*nz ) * Math.sqrt( bx*bx + by*by + bz*bz ) ) );
				}

				if ( Math.abs( ev ) > 0.01 )
				{
					// set attenuated intensity
					rOut.get().set( (float)Math.toDegrees( dotP ) );
				}

				// TODO: update ray position

				final double v = rIn.get().get();

				rIn.bck( 1 );
				rOut.bck( 1 );
			}
		}

		return img;
	}

	public static Img< FloatType > simulate()
	{
		return simulate( false, rnd );
	}

	public static Img< FloatType > simulate( final boolean halfPixelOffset, final Random rnd )
	{
		// rendering in a higher resolution and then downsampling it makes
		// it much smoother and somewhat more realistic, otherwise there
		// are artifacts from linear interpolation and hard borders
		int scale = 2;
		int size = 289;
		
		if ( scale == 2 )
			size++; // one pixel is lost when downsampling
		
		// open with ImgOpener using an ImagePlusImg
        Img< FloatType > img = new ArrayImgFactory< FloatType >().create( new long[] { size*scale, size*scale, size*scale }, new FloatType() );

        // draw a small sphere for every pixel of a larger sphere
        drawSpheres( img, 0, 1, scale, halfPixelOffset, rnd );
        
        if ( scale == 2 )
        	img = downSample2x( img );
        
		return img;
	}
	
	public static Img< FloatType > downSample2x( final RandomAccessibleInterval< FloatType > randomAccessible )
	{
		final long[] dim = new long[ randomAccessible.numDimensions() ];
		
		for ( int d = 0; d < dim.length; ++d )
			dim[ d ] = randomAccessible.dimension( d ) / 2 - 1;
		
		final Img< FloatType > img = new ArrayImgFactory< FloatType >().create( dim, new FloatType() );
		
		final NLinearInterpolatorFactory< FloatType > factory = new NLinearInterpolatorFactory< FloatType >();
		final RealRandomAccessible< FloatType > interpolant = Views.interpolate( Views.extendMirrorSingle( randomAccessible ), factory );
		final RealRandomAccess< FloatType > realRandomAccess = interpolant.realRandomAccess();
		
		final Cursor< FloatType > c = img.localizingCursor();
		final int[] l = new int[ img.numDimensions() ];
		final double[] lf = new double[ img.numDimensions() ];
		
		while ( c.hasNext() )
		{
			c.fwd();
			c.localize( l );

			for ( int d = 0; d < dim.length; ++d )
				lf[ d ] = l[ d ]*2.0 + 0.5;

			realRandomAccess.setPosition( lf );
			c.get().set( realRandomAccess.get() );
		}
		
		return img;
	}
	
	/**
     * Draws a sphere that contains lots of small spheres into the center of the interval
     *
     * @param randomAccessible - the image data to write to
     * @param minValue - the minimal intensity of one of the small spheres
     * @param maxValue - the maximal intensity of one of the small spheres
     * @param scale - the scale
     * @param rnd - the random number generator
     * @param <T> - the type
     */
    public static < T extends RealType< T > > void drawSpheres(
            final RandomAccessibleInterval< T > randomAccessible,
            final double minValue, final double maxValue, final int scale,
            final boolean halfPixelOffset,
            final Random rnd )
    {
            // the number of dimensions
            int numDimensions = randomAccessible.numDimensions();

            // define the center and radius
            Point center = new Point( randomAccessible.numDimensions() );
            long minSize = randomAccessible.dimension( 0 );

            for ( int d = 0; d < numDimensions; ++d )
            {
                    long size = randomAccessible.dimension( d );

                    center.setPosition( size / 2 , d );
                    minSize = Math.min( minSize, size );
            }

            // define the maximal radius of the small spheres
            int maxRadius = 10 * scale;

            // compute the radius of the large sphere so that we do not draw
            // outside of the defined interval
            long radiusLargeSphere = minSize / 2 - 47*scale - 1;

            // instantiate a random number generator
            //Random rnd = new Random( System.currentTimeMillis() );

            // define a hypersphere (n-dimensional sphere)
            HyperSphere< T > hyperSphere =
                    new HyperSphere<T>( randomAccessible, center, radiusLargeSphere );

            // create a cursor on the hypersphere
            HyperSphereCursor< T > cursor = hyperSphere.cursor();

            int size = (int)hyperSphere.size();

            IJ.showProgress( 0.0 );

            int i = 0;

            while ( cursor.hasNext() )
            {
                    cursor.fwd();

                    // the random radius of the current small hypersphere
                    int radius = rnd.nextInt( maxRadius ) + 1;

                    // instantiate a small hypersphere at the location of the current pixel
                    // in the large hypersphere
                    final HyperSphere< T > smallSphere;
                    
                    if ( halfPixelOffset )
                    {
                    	// shifting by one pixel in xy means half a pixel after downsampling
                    	final long[] tmp = new long[ randomAccessible.numDimensions() ];
                    	cursor.localize( tmp );
                    	for ( int d = 0; d < tmp.length - 1; ++d )
                    		tmp[ d ] += 1;
                    	smallSphere = new HyperSphere< T >( randomAccessible, new Point( tmp ), radius );
                    }
                    else
                    {
                    	smallSphere = new HyperSphere< T >( randomAccessible, cursor, radius );
                    }

                    // define the random intensity for this small sphere
                    double randomValue = rnd.nextDouble();

                    // take only every 4^dimension'th pixel by chance so that it is not too crowded
                    if ( Math.round( randomValue * 10000 ) % Util.pow( 7*scale, numDimensions ) == 0 )
                    {
                            // scale to right range
                            randomValue = rnd.nextDouble() * ( maxValue - minValue ) + minValue;

                            // set the value to all pixels in the small sphere if the intensity is
                            // brighter than the existing one
                            for ( final T value : smallSphere )
                                    value.setReal( Math.max( randomValue, value.getRealDouble() ) );
                    }

                    IJ.showProgress( ++i, size );
            }
    }

	public static void main( String[] args )
	{
		final String dir = "src/main/resources/";
		final ExecutorService service = Executors.newFixedThreadPool( Runtime.getRuntime().availableProcessors() );

		new ImageJ();

		final float poissonSNR = 25f;
		final int lightsheetSpacing = 3;
		final float attenuation = 0.01f;
		
		// six angles
		//final int angleIncrement = 60;
		//final float osem = 3.0f;
		
		// seven angles
		final int angleIncrement = 52;
		final float osem = 3.0f;

		// eight angles
		//final int angleIncrement = 45;
		//final float osem = 4.0f;

		// so that everything is rotated at least once
		final int angleOffset = 15;
		
		// artificially rendered object based on which everything is computed,
		// including the ground-truth image, which is rotated once by the angle offset
		// so that it is realistic
		System.out.println( new Date( System.currentTimeMillis() ) + ": rendering basis for ground truth" );
		final Img<FloatType> rendered = simulate();

		//System.out.println( new Date( System.currentTimeMillis() ) + ": computing ground truth" );
		//final Img<FloatType> obj = rotateAroundAxis( rendered, 0, angleOffset );

		//Tools.save( rendered, dir + "rendered.tif" );
		//Tools.save( obj, dir + "groundtruth.tif" );

		final ArrayList< Img< FloatType > > weights = new ArrayList< Img< FloatType > >();
		
		//ImageJFunctions.show( rendered ).setTitle( "rendered" );
		//ImageJFunctions.show( obj ).setTitle( "groundtruth" );
		
		for ( int angle = 0; angle < 360; angle += angleIncrement )
		{
			System.out.println( new Date( System.currentTimeMillis() ) + ": rotating angle " + angle );
			Img<FloatType> rot = rotateAroundAxis( rendered, 0, angle + angleOffset );

			System.out.println( new Date( System.currentTimeMillis() ) + ": refracting angle " + angle );

			drawSimpleImage( rot );
			Pair< Img< FloatType >, Img< FloatType > > eigen = Hessian.largestEigenVector( rot );
			Img<FloatType> refr = refract3d( rot );
			
			ImageJFunctions.show( rot ).setDisplayRange( 0, 1 );
			ImageJFunctions.show( eigen.getA() ).setDisplayRange( -1, 1 );
			ImageJFunctions.show( eigen.getB() ).setDisplayRange( -1, 1 );
			ImageJFunctions.show( refr ).setDisplayRange( 0, 90 );
			SimpleMultiThreading.threadHaltUnClean();

			System.out.println( new Date( System.currentTimeMillis() ) + ": attenuation angle " + angle );
			Img<FloatType> att = attenuate3d( rot, attenuation );

			System.out.println( new Date( System.currentTimeMillis() ) + ": weights angle " + angle );
			Img<FloatType> w = computeWeightImage( rot, attenuation );
			
			System.out.println( new Date( System.currentTimeMillis() ) + ": convolving angle " + angle );
			Img< FloatType > psf = Tools.open( dir + "Angle" + angle + ".tif", true );
			Img<FloatType> con = convolve( att, psf, service );
			
			Tools.adjustImage( con, minValue, avgIntensity );
			
			System.out.println( new Date( System.currentTimeMillis() ) + ": extracting slices angle " + angle );
			Img<FloatType> acq = extractSlices( con, lightsheetSpacing, poissonSNR );
			
			System.out.println( new Date( System.currentTimeMillis() ) + ": make isotropic angle " + angle );
			Img<FloatType> iso = makeIsotropic( acq, lightsheetSpacing );
			
			System.out.println( new Date( System.currentTimeMillis() ) + ": rotate back view, weights and psf angle " + angle );
			Img<FloatType> view = rotateAroundAxis( iso, 0, -angle );	
			Img<FloatType> viewWeights = rotateAroundAxis( w, 0, -angle );			
			Img<FloatType> viewPSF = rotateAroundAxis( psf, 0, -angle );
			
			//ImageJFunctions.show( view ).setTitle( "aligned_view_" + angle );
			//ImageJFunctions.show( viewPSF ).setTitle( "aligned_view_psf_" + angle );
			
			Tools.save( rot, dir + "rot_view_" + angle + ".tif" );
			Tools.save( att, dir + "att_view_" + angle + ".tif" );
			Tools.save( con, dir + "con_view_" + angle + ".tif" );
			Tools.save( acq, dir + "acq_view_" + angle + ".tif" );
			Tools.save( iso, dir + "iso_view_" + angle + ".tif" );
			Tools.save( view, dir + "aligned_view_" + angle + ".tif" );
			Tools.save( viewPSF, dir + "aligned_view_psf_" + angle + ".tif" );
			
			weights.add( viewWeights );
			//ImageJFunctions.show( att ).setTitle( "attenuated" );
			//ImageJFunctions.show( con ).setTitle( "convolved" );
			//ImageJFunctions.show( acq ).setTitle( "acquisition" );
			//ImageJFunctions.show( iso ).setTitle( "isotropic" );
			
			//ImageJFunctions.show( viewWeights ).duplicate().show();
		}
		
		// norm sum weights to 3
		final ArrayList< Cursor< FloatType > > cursors = new ArrayList< Cursor<FloatType > >(); 

		for ( final Img< FloatType > w : weights )
			cursors.add( w.cursor() );
		
		final Cursor< FloatType > c = cursors.get( 0 );
		
		while ( c.hasNext() )
		{
			float sum = 0;
			
			for ( final Cursor< FloatType > cursor : cursors )
				sum += cursor.next().get();
			
			if ( sum == 0 )
			{
				for ( final Cursor< FloatType > cursor : cursors )
					cursor.get().set( 0 );
			}
			else
			{
				for ( final Cursor< FloatType > cursor : cursors )
					cursor.get().set( Math.min( 1, osem * (cursor.get().get() / (float)sum) ) );
			}
		}
		
		for ( int i = 0; i < weights.size(); ++i )
		{
			final int angle = i * angleIncrement;
			Tools.save( weights.get( i ), dir + "aligned_view_weights" + angle + ".tif" );	
		}
		
		final Img< FloatType > sumWeights = weights.get( 0 ).factory().create( weights.get( 0 ), weights.get( 0 ).firstElement() );
		
		for ( final Img< FloatType > w : weights )
		{
			final Cursor< FloatType > s = sumWeights.cursor();
			final Cursor< FloatType > v = w.cursor();
			
			while( s.hasNext() )
			{
				s.fwd();
				v.fwd();
				s.get().add( v.get() );
			}
		}
		
		Tools.save( sumWeights, dir + "sum_weights.tif" );
		System.out.println( "done" );
	}
}
