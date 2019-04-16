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

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.IJ;
import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.region.hypersphere.HyperSphereCursor;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.multithreading.SimpleMultiThreading;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.Util;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;
import net.preibisch.simulation.cluster.RunJob;
import net.preibisch.simulation.raytracing.Lightsheet;
import net.preibisch.simulation.raytracing.Raytrace;

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

	private static final boolean inside( final double[] rayPosition, final Interval interval )
	{
		for ( int d = 0; d < rayPosition.length; ++d )
			if ( rayPosition[ d ] < interval.min( d ) || rayPosition[ d ] > interval.max( d ) )
				return false;

		return true;
	}

	public static Img< FloatType > projectToCamera(
			//final RandomAccessibleInterval< FloatType > imgIn,
			final RandomAccessibleInterval< FloatType > imgRi,
			final RandomAccessibleInterval< FloatType > refr, final double ri, final int currentzPlane )
	{
		// the refracted image
		final Img< FloatType > proj = new ArrayImgFactory< FloatType >( new FloatType() ).create( new long[] { refr.dimension( 0 ), refr.dimension( 1 ) } );

		final Img< FloatType > tmp = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgRi );
		final RandomAccess< FloatType > rt = Views.extendZero( tmp ).randomAccess();

		//final RealRandomAccess< FloatType > rrIm = Views.interpolate( Views.extendMirrorSingle( imgIn ), new NLinearInterpolatorFactory<>() ).realRandomAccess();
		final RealRandomAccess< FloatType > rrRi = Views.interpolate( Views.extendMirrorSingle( imgRi ), new NLinearInterpolatorFactory<>() ).realRandomAccess();
		final RealRandomAccess< FloatType > rrRef = Views.interpolate( Views.extendMirrorSingle( refr ), new NLinearInterpolatorFactory<>() ).realRandomAccess();

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		// the incoming direction of the light
		final double[] rayVector = new double[ 3 ];

		// the outgoing, refracted direction of the light
		final double[] refractedRay = new double[ 3 ];

		// the current position of the ray
		final double[] rayPosition = new double[ 3 ];

		final double nA = 1.00; // air (intensity == 0)
		final double nB = 1.01; //ri // water (intensity == 1)

		final Cursor< FloatType > c = proj.localizingCursor();

		final double sigma = 4.0;
		final double two_sq_sigma = 2 * sigma * sigma;

		while ( c.hasNext() )
		{
			c.fwd();

			final boolean debug = false;

			//if ( c.getIntPosition( 0 ) == 78 && c.getIntPosition( 1 ) == 75 )
			//	debug = true;

			double avgValue = 0;

			for ( int i = 0; i < 500; ++i )
			{
				rayPosition[ 0 ] = c.getIntPosition( 0 ) + (rnd.nextDouble() - 0.5);
				rayPosition[ 1 ] = c.getIntPosition( 1 ) + (rnd.nextDouble() - 0.5);
				rayPosition[ 2 ] = 1;
	
				rayVector[ 0 ] = 0;//(rnd.nextDouble() - 0.5)/100;
				rayVector[ 1 ] = 0;//(rnd.nextDouble() - 0.5)/100;
				rayVector[ 2 ] = 1;
	
				double signal = 0;
				long maxMoves = imgRi.dimension( 2 );
				int moves = 0;

				//if ( debug )
					//System.out.println( "\n" + Util.printCoordinates( rayPosition ) );

				while ( inside( rayPosition, refr ) && moves < maxMoves )
				{
					++moves;
					rrRi.setPosition( rayPosition );
	
					// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)
					Hessian.computeHessianMatrix3D( rrRi, matrix );
	
					double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );
	
					//if ( debug )
						//System.out.println( "ev: " + ev );

					if ( Math.abs( ev ) > 0.01 )
					{
						// compute refractive index change
						rrRi.setPosition( rayPosition );
	
						rrRi.move( -rayVector[ 0 ], 0 );
						rrRi.move( -rayVector[ 1 ], 1 );
						rrRi.move( -rayVector[ 2 ], 2 );
	
						// intensity at the origin of the ray
						final double i0 = rrRi.get().get();
	
						rrRi.move( 2*rayVector[ 0 ], 0 );
						rrRi.move( 2*rayVector[ 1 ], 1 );
						rrRi.move( 2*rayVector[ 2 ], 2 );
	
						// intensity at the projected location of the ray
						final double i1 = rrRi.get().get();
	
						final double n0 = ( nB - nA ) * i0 + nA;
						final double n1 = ( nB - nA ) * i1 + nA;

						//if ( debug )
						//	System.out.println( n0 + " " + n1 );

						final double thetaI = Raytrace.incidentAngle( rayVector, eigenVector );
						final double thetaT = Raytrace.refract( rayVector, eigenVector, n0, n1, thetaI, refractedRay );
	
						// total reflection
						if ( Double.isNaN( thetaT ) )
						{
							if ( debug )
								System.out.println( "reflect: " + thetaI + " " + thetaT );

							//Raytrace.reflect( rayVector, eigenVector, refractedRay );
							refractedRay[ 0 ] = rayVector[ 0 ];
							refractedRay[ 1 ] = rayVector[ 1 ];
							refractedRay[ 2 ] = rayVector[ 2 ];
						}

						Raytrace.norm( refractedRay );
	
						// update the ray vector
						rayVector[ 0 ] = refractedRay[ 0 ];
						rayVector[ 1 ] = refractedRay[ 1 ];
						rayVector[ 2 ] = refractedRay[ 2 ];

						if ( debug )
							System.out.println( Util.printCoordinates( rayVector ) );
					}
	
					// place a gaussian sphere
					rrRef.setPosition( rayPosition );
					final double zOffset = Math.abs( rayPosition[ 2 ] - currentzPlane );
					/*
					 * 0 > 1
					 * 0.5 > 1.5
					 * 1 > 2
					 */
					//signal += ( rrRef.get().get() / Math.pow( zOffset + 1.0, 1.00 / 1.75 ));
					signal += rrRef.get().get() * getGaussValue( zOffset, two_sq_sigma );
	
					rayPosition[ 0 ] += rayVector[ 0 ];
					rayPosition[ 1 ] += rayVector[ 1 ];
					rayPosition[ 2 ] += rayVector[ 2 ];

					if ( debug )
					{
						rt.setPosition( Math.round( rayPosition[ 0 ] ), 0 );
						rt.setPosition( Math.round( rayPosition[ 1 ] ), 1 );
						rt.setPosition( Math.round( rayPosition[ 2 ] ), 2 );
						rt.get().set( rt.get().get() + 1.0f );

						//System.out.println( "signal: " + signal + " @ " + Util.printCoordinates( rayPosition ) );
					}
				}

				avgValue += signal;

			}

			c.get().set( (float)(avgValue/10.0) );

			if ( debug )
				ImageJFunctions.show( tmp );

		}

		return proj;
	}

	private static final double getGaussValue( final double distance, final double two_sq_sigma )
	{
		return Math.exp( -(distance * distance) / two_sq_sigma );
	}

	public static VolumeInjection refract3d(
			final RandomAccessibleInterval< FloatType > imgIn,
			final RandomAccessibleInterval< FloatType > imgRi,
			final boolean illum, final int z, final double lsMiddle, final double lsEdge, final double ri )
	{
		// the refracted image
		final Img< FloatType > img = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );
		final Img< FloatType > weight = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );

		// for drawing individual rays ...
		//final RandomAccess< FloatType > rOut = img.randomAccess();
		//final long[] intPos = new long[ 3 ];

		final RealRandomAccess< FloatType > rrIm = Views.interpolate( Views.extendMirrorSingle( imgIn ), new NLinearInterpolatorFactory<>() ).realRandomAccess();
		final RealRandomAccess< FloatType > rrRi = Views.interpolate( Views.extendMirrorSingle( imgRi ), new NLinearInterpolatorFactory<>() ).realRandomAccess();

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		final double[] sigma = new double[] { 0.5, 0.5, 0.5 };
		final VolumeInjection inject = new VolumeInjection( img, weight, sigma );

		// the incoming direction of the light
		final double[] rayVector = new double[ 3 ];

		// the outgoing, refracted direction of the light
		final double[] refractedRay = new double[ 3 ];

		// the current position of the ray
		final double[] rayPosition = new double[ 3 ];

		final double nA = 1.00; // air (intensity == 0)
		final double nB = ri; // water (intensity == 1)

		System.out.println( "middle: " + lsMiddle );
		System.out.println( "lsEdge: " + lsEdge );
		final Lightsheet ls = new Lightsheet( img.dimension( 0 ) / 2.0, lsMiddle, img.dimension( 0 ), lsEdge );

		System.out.println( "starting..." );
		long time = System.currentTimeMillis();

		// for each point on the xz plane at y=0
		final int numRays = 200000;
		final long maxMoves = img.dimension( 2 );
		final Random rnd = new Random( 2423 );

		for ( int i = 0; i < numRays; ++i )
		{
			rayPosition[ 0 ] = rnd.nextDouble() * imgIn.max( 0 ); //x;//c.getIntPosition( 0 );
			rayPosition[ 1 ] = illum ? (int)imgIn.dimension( 1 ) - 1 : 0;

			final double lightsheetthickness = ls.predict( rayPosition[ 0 ] );
			rayPosition[ 2 ] = z + (rnd.nextDouble()*lightsheetthickness)-lightsheetthickness/2.0;//c.getIntPosition( 1 );

			rayVector[ 0 ] = (rnd.nextDouble() - 0.5)/5;
			rayVector[ 1 ] = illum ? -1 : 1;
			rayVector[ 2 ] = 0;

			Raytrace.norm( rayVector );

			int moves = 0;

			while ( inside( rayPosition, imgIn ) && moves < maxMoves )
			{
				++moves;

				rrIm.setPosition( rayPosition );
				rrRi.setPosition( rayPosition );
				final float valueIm = rrIm.get().get();

				// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)
				Hessian.computeHessianMatrix3D( rrRi, matrix );

				double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

				if ( Math.abs( ev ) > 0.01 )
				{
					// compute refractive index change
					rrRi.setPosition( rayPosition );

					rrRi.move( -rayVector[ 0 ], 0 );
					rrRi.move( -rayVector[ 1 ], 1 );
					rrRi.move( -rayVector[ 2 ], 2 );

					// intensity at the origin of the ray
					final double i0 = rrRi.get().get();

					rrRi.move( 2*rayVector[ 0 ], 0 );
					rrRi.move( 2*rayVector[ 1 ], 1 );
					rrRi.move( 2*rayVector[ 2 ], 2 );

					// intensity at the projected location of the ray
					final double i1 = rrRi.get().get();

					final double n0 = ( nB - nA ) * i0 + nA;
					final double n1 = ( nB - nA ) * i1 + nA;

					final double thetaI = Raytrace.incidentAngle( rayVector, eigenVector );
					final double thetaT = Raytrace.refract( rayVector, eigenVector, n0, n1, thetaI, refractedRay );

					// total reflection
					if ( Double.isNaN( thetaT ) )
					{
						refractedRay[ 0 ] = rayVector[ 0 ];
						refractedRay[ 1 ] = rayVector[ 1 ];
						refractedRay[ 2 ] = rayVector[ 2 ];

						//continue; // that's an expensive getting stuck
						//Raytrace.reflect( rayVector, eigenVector, refractedRay );
					}

					Raytrace.norm( refractedRay );

					// update the ray vector
					rayVector[ 0 ] = refractedRay[ 0 ];
					rayVector[ 1 ] = refractedRay[ 1 ];
					rayVector[ 2 ] = refractedRay[ 2 ];
				}

				// place a gaussian sphere
				inject.addNormalizedGaussian( valueIm, rayPosition );

				/*
				// for drawing individual rays ...
				intPos[ 0 ] = Math.round( rayPosition[ 0 ] );
				intPos[ 1 ] = Math.round( rayPosition[ 1 ] );
				intPos[ 2 ] = Math.round( rayPosition[ 2 ] );
				rOut.setPosition( intPos );
				rOut.get().set( value );
				*/

				rayPosition[ 0 ] += rayVector[ 0 ];
				rayPosition[ 1 ] += rayVector[ 1 ];
				rayPosition[ 2 ] += rayVector[ 2 ];
			}
		}

		System.out.println( " ... " + ( System.currentTimeMillis() - time ) );

		return inject;
	}

	public static Pair< Img< FloatType >, Img< FloatType > > simulate( final String dir)
	{
		return simulate( rnd, dir );
	}

	public static Pair< Img< FloatType >, Img< FloatType > > simulate( final Random rnd, final String dir )
	{
		// rendering in a higher resolution and then downsampling it makes
		// it much smoother and somewhat more realistic, otherwise there
		// are artifacts from linear interpolation and hard borders
		int scale = 2;
		int size = 289;
		
		if ( scale == 2 )
			size++; // one pixel is lost when downsampling
		
		// open with ImgOpener using an ImagePlusImg
        Img< FloatType > img = new ArrayImgFactory< FloatType >( new FloatType() ).create( new long[] { size*scale, size*scale, size*scale } );
        System.out.println( "Loading " + dir + "block4.tif" );
        Img< FloatType > ri = Tools.open( dir + "block4.tif", new ArrayImgFactory< FloatType >( new FloatType() ) );

        System.out.println( "Adding noise" );
        for ( final FloatType t : ri )
        	t.setReal( Math.max( 0, t.get() + (rnd.nextDouble() - 0.5)/10 ) );
        //Gauss3.gauss( 3, Views.extendMirrorSingle( img ), img );

        // draw a small sphere for every pixel of a larger sphere
        System.out.println( "Adding spheres" );
        multiSpheres( img, ri, scale, rnd );
        
        if ( scale == 2 )
        {
        	img = downSample2x( img );
        	ri = downSample2x( ri );
        }
        
		return new ValuePair<>( img, ri );
	}
	
	public static Img< FloatType > downSample2x( final RandomAccessibleInterval< FloatType > randomAccessible )
	{
		final long[] dim = new long[ randomAccessible.numDimensions() ];
		
		for ( int d = 0; d < dim.length; ++d )
			dim[ d ] = randomAccessible.dimension( d ) / 2 - 1;
		
		final Img< FloatType > img = new ArrayImgFactory< FloatType >( new FloatType() ).create( dim );
		
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

	public static < T extends RealType< T > > void multiSpheres(
			final RandomAccessibleInterval< T > image,
			final RandomAccessibleInterval< T > ri,
			final int scale,
			final Random rnd )
	{
			// the number of dimensions
			final int numDimensions = image.numDimensions();

			// define the center and radius
			final Point center1 = new Point( image.numDimensions() );
			final Point center2 = new Point( image.numDimensions() );
			final Point center3 = new Point( image.numDimensions() );

			long minSize = image.dimension( 0 );

			for ( int d = 0; d < numDimensions; ++d )
			{
				long size = image.dimension( d );

				center1.setPosition( size / 2 , d );
				center2.setPosition( size / 2 , d );
				center3.setPosition( size / 2 , d );
				minSize = Math.min( minSize, size );
			}

			// define the maximal radius of the small spheres
			int maxRadius = 10 * scale;

			// compute the radius of the large sphere so that we do not draw
			// outside of the defined interval
			long radiusLargeSphere1 = minSize / 2 - 47*scale - 1 /*new*/;// -45;

			final ArrayList< Pair< Pair< Point, Long >, double[] > > centers = new ArrayList<>();

			// defines a small sphere (location, radius), [minValueIm, maxValueIm, minValueRi, maxValueRi]
			centers.add( new ValuePair<>( new ValuePair<>( center2, radiusLargeSphere1 ), new double[] { 0.5, 1.0, 1.0, 1.1 } ) );

			for ( final Pair< Pair< Point, Long >, double[] > center : centers )
			{
				final double minValueIm = center.getB()[ 0 ];
				final double maxValueIm = center.getB()[ 1 ];

				final double minValueRi = center.getB()[ 2 ];
				final double maxValueRi = center.getB()[ 3 ];

				// define a hypersphere (n-dimensional sphere)
				final HyperSphere< T > hyperSphereIm = new HyperSphere<T>( image, center.getA().getA(), center.getA().getB() );
				final HyperSphere< T > hyperSphereRi = new HyperSphere<T>( ri, center.getA().getA(), center.getA().getB() );

				// create a cursor on the hypersphere
				final HyperSphereCursor< T > cursorIm = hyperSphereIm.cursor();
				final HyperSphereCursor< T > cursorRi = hyperSphereRi.cursor();
				
				final int size = (int)hyperSphereIm.size();
				
				IJ.showProgress( 0.0 );
				
				int i = 0;
				
				while ( cursorIm.hasNext() )
				{
					cursorIm.fwd();
					cursorRi.fwd();

					if ( maxValueIm == minValueIm )
							cursorIm.get().setReal( minValueIm );

					if ( maxValueRi == minValueRi )
							cursorRi.get().setReal( minValueRi );

					// the random radius of the current small hypersphere
					int radius = Math.max( rnd.nextInt( maxRadius ) + 1, maxRadius - 1 );

					// instantiate a small hypersphere at the location of the current pixel
					// in the large hypersphere
					final HyperSphere< T > smallSphereIm = new HyperSphere< T >( image, cursorIm, radius );
					final HyperSphere< T > smallSphereRi = new HyperSphere< T >( ri, cursorRi, radius );
					
					// define the random intensity for this small sphere
					double randomValue = rnd.nextDouble();
					
					// take only every 4^dimension'th pixel by chance so that it is not too crowded
					if ( ( randomValue * 100000 ) < 1 )
					{
						// scale to right range
						randomValue = rnd.nextDouble();
						
						double randomValueRi = randomValue * ( maxValueRi - minValueRi ) + minValueRi;
						double randomValueIm = randomValue * ( maxValueIm - minValueIm ) + minValueIm;
						
						// set the value to all pixels in the small sphere if the intensity is
						// brighter than the existing one
						if ( maxValueIm != minValueIm )
							for ( final T value : smallSphereIm )
								value.setReal( Math.max( randomValueIm, value.getRealDouble() ) );
						
						if ( maxValueRi != minValueRi )
						{
							for ( final T value : smallSphereRi )
							{
								if ( value.getRealDouble() == 5.0 )
									value.setReal( randomValueRi );
								else
									value.setReal( Math.max( randomValueRi, value.getRealDouble() ) );
							}
						}
					}

					IJ.showProgress( ++i, size );
				}
			}
	}

	public static void simulate( final boolean illum, final double lsMiddle, final double lsEdge, final double ri, final String dir, final ExecutorService service, final int z )
	{
		final ArrayList< Integer > zPlanes = new ArrayList<>();
		zPlanes.add( z );
		simulate( illum, lsMiddle, lsEdge, ri, dir, service, zPlanes );
	}

	public static void simulate( final boolean illum, final double lsMiddle, final double lsEdge, final double ri, final String dir, final ExecutorService service )
	{
		simulate( illum, lsMiddle, lsEdge, ri, dir, service, null );
	}

	public static void simulate( final boolean illum, final double lsMiddle, final double lsEdge, final double ri, final String dir, final ExecutorService service, ArrayList< Integer > zPlanes )
	{
		// artificially rendered object based on which everything is computed,
		// including the ground-truth image, which is rotated once by the angle offset
		// so that it is realistic
		System.out.println( new Date( System.currentTimeMillis() ) + ": rendering basis for ground truth" );
		final Pair< Img<FloatType>, Img<FloatType> > rendered = simulate( new Random( 464232194 ), dir );

		if ( zPlanes == null )
		{
			zPlanes = new ArrayList<>();
			for ( int z = 0; z < rendered.getA().dimension( 2 ); ++z )
				zPlanes.add( z );
		}

		//Tools.save( rendered, dir + "rendered.tif" );

		RunJob.display( rendered.getA(), "rendered" );
		RunJob.display( rendered.getB(), "RI_rendered" );
		//SimpleMultiThreading.threadHaltUnClean();

		Img<FloatType> rotIm = rendered.getA();
		Img<FloatType> rotRi = rendered.getB();

		System.out.println( new Date( System.currentTimeMillis() ) + ": refracting ... " );

		for ( final int z : zPlanes )
		{
			String tag = illum + "_" + ri + "_" + z;
			System.out.println( tag );
			System.out.println( new File(  dir + "refr_img_" + tag + ".tif" ).getAbsolutePath() );

			if ( !new File(  dir + "refr_img_" + tag + ".tif" ).exists() )
			{
				VolumeInjection simulated = refract3d( rotIm, rotRi, illum, z, lsMiddle, lsEdge, ri );
				Tools.save( simulated.getImage(), dir + "refr_img_" + tag + ".tif" );
				Tools.save( simulated.getWeight(), dir + "refr_weight_" + tag + ".tif" );
			}

			//SimpleMultiThreading.threadHaltUnClean();

			//RunJob.display( refr.getImage(), "img" ).show();
			//RunJob.display( refr.getWeight(), "weight" ).show();
			//RunJob.display( refr.normalize(), "normed" ).show();

			Img<FloatType> refr = Tools.open( dir + "refr_img_" + tag + ".tif", new ArrayImgFactory<>() );
			//Img<FloatType> weight = Tools.open( dir + "refr_weight_" + tag + ".tif", new ArrayImgFactory<>() );
			//Img<FloatType> refr = VolumeInjection.normalize( image, weight );
			RunJob.display( refr, "refr" );
			//RunJob.display( weight, "weight" );
			//RunJob.display( refr, "norm" );

			Img< FloatType> proj = projectToCamera( rotRi, refr, ri, z );
			//RunJob.display( proj, "proj" );
			Tools.save( proj, dir + "proj_" + tag + ".tif" );
		}

		//ImageJFunctions.show( rot ).setDisplayRange( 0, 1 );
		//ImageJFunctions.show( eigen.getA() ).setDisplayRange( -1, 1 );
		//ImageJFunctions.show( eigen.getB() ).setDisplayRange( -1, 1 );
		//ImageJFunctions.show( refr ).setDisplayRange( 0, 1 );
		//if ( RunJob.isCluster )
		//	System.exit( 0 );
		//else
		//	SimpleMultiThreading.threadHaltUnClean();


		System.out.println( "done" );
	}

	public static void main( String[] args )
	{
		final double sigma = 4;
		final double two_sq_sigma = 2 * sigma * sigma;

		for ( double d = 0; d < 20; d += 1 )
			System.out.println( d + ": " + getGaussValue( d, two_sq_sigma ) );

		//System.exit( 0 );

		final double[] i = new double[] { 0, -1, 0 };
		final double[] n = new double[] { 0.7071067811865475,-0.7071067811865475,-0.0 }; // 45
		//final double[] n = new double[] { 0.642824346533225,-0.766013615743305,-0.0 }; // 40.00274776305653

		// refraction
		final double[] t = new double[ 3 ];

		// relection
		final double[] r = new double[ 3 ];

		// 1.0 > 1.2: 45 > 36.10420471349619
		// 1.0 > 1.1: 45 > 40.00274776305653
		// 1.1 > 1.2: 40.00274776305653 > 36.10420471349619
		final double n0 = 1.1;
		final double n1 = 1.2;

		final double thetaI = Raytrace.incidentAngle( i, n );
		final double thetaT = Raytrace.refract( i, n, n0, n1, thetaI, t );
		Raytrace.reflect( i, n, r );

		System.out.println( "Refraction: " );
		System.out.println( Math.toDegrees( thetaI ) + " >> " + Math.toDegrees( thetaT ) );
		System.out.println( "i: " + Util.printCoordinates( i ) );
		System.out.println( "n: " + Util.printCoordinates( n ) );
		System.out.println( "t: " + Util.printCoordinates( t ) );

		System.out.println( "Reflection: " );
		System.out.println( "t: " + Util.printCoordinates( r ) );

		final String dir = "render/";
		final ExecutorService service = Executors.newFixedThreadPool( Runtime.getRuntime().availableProcessors() );

		new ImageJ();

		//ImageJFunctions.show( Tools.openBF( dir + "block3.tif", new ArrayImgFactory< FloatType >() ) );
		//SimpleMultiThreading.threadHaltUnClean();

		final int z = 176;
		final double lsMiddle = 1.0;
		final double lsEdge = 3.0;
		final double ri = 1.1;
		final boolean illum = false; // true > bottom, false > top

		simulate( illum, lsMiddle, lsEdge, ri, dir, service, z );
		//simulate( !illum, lsMiddle, lsEdge, ri, dir, service, z );
	}
}
