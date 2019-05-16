package net.preibisch.simulation.raytracing;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.TextRoi;
import ij.gui.Toolbar;
import ij.process.FloatProcessor;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import net.preibisch.simulation.Hessian;
import net.preibisch.simulation.SimulateMultiViewAberrations;
import net.preibisch.simulation.VolumeInjection;

public class RayTracingTest
{
	public static Img< FloatType > hessian( final Img< FloatType > img, final int z )
	{
		final Img< FloatType > hessian = img.copy();

		final Cursor< FloatType > cursor = img.localizingCursor();
		final RandomAccess< FloatType > ra = Views.extendMirrorSingle( img ).randomAccess();
		final RandomAccess< FloatType > raH = hessian.randomAccess();

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		while( cursor.hasNext() )
		{
			cursor.fwd();
			if ( cursor.getIntPosition( 2 ) != z )
				continue;

			ra.setPosition( cursor );
			Hessian.computeHessianMatrix3D( ra, matrix );
			final double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

			if ( cursor.getIntPosition( 0 ) == 85 && cursor.getIntPosition( 1 ) == 85 )
				System.out.println( Util.printCoordinates( eigenVector ) );

			raH.setPosition( cursor );
			raH.get().set( (float) ev);
		}

		return hessian;
	}

	public static VolumeInjection drawMultipleRays( final RandomAccessibleInterval< FloatType > imgIn, final int z, final int x, final int numRays )
	{
		// the refracted image
		final Img< FloatType > imgOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );
		final Img< FloatType > weightOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );

		final double[] sigma = new double[] { 0.5, 0.5, 0.5 };
		final VolumeInjection inject = new VolumeInjection( imgOut, weightOut, sigma );

		// the incoming direction of the light
		final double[] rayVector = new double[ 3 ];

		// the outgoing, refracted direction of the light
		final double[] refractedRay = new double[ 3 ];

		// the current position of the ray
		final double[] rayPosition = new double[ imgIn.numDimensions() ];

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		//final RandomAccess< FloatType > ra = img.randomAccess();
		final RealRandomAccess< FloatType > rr = Views.interpolate( Views.extendMirrorSingle( imgIn ), new NLinearInterpolatorFactory<>() ).realRandomAccess();

		final Random rnd = new Random( 234345 );
		final double valueIm = 1.0;

		for ( int n = 0; n < numRays; ++n )
		{
			rayVector[ 0 ] = 0;
			rayVector[ 1 ] = -1;
			rayVector[ 2 ] = 0;

			Raytrace.norm( rayVector );

			rayPosition[ 0 ] = x + rnd.nextDouble() - 0.5;
			rayPosition[ 1 ] = imgIn.dimension( 1 ) - 1;
			rayPosition[ 2 ] = z + rnd.nextDouble() - 0.5;

			while ( SimulateMultiViewAberrations.inside( rayPosition, imgIn ) )
			{
				rr.setPosition( rayPosition );

				// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)
				Hessian.computeHessianMatrix3D( rr, matrix );

				double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

				if ( Math.abs( ev ) > 0.01 )
				{
					// compute refractive index change
					rr.setPosition( rayPosition );

					rr.move( -rayVector[ 0 ], 0 );
					rr.move( -rayVector[ 1 ], 1 );
					rr.move( -rayVector[ 2 ], 2 );

					// refractive index at the origin of the ray
					final double n0 = rr.get().get();

					rr.move( 2*rayVector[ 0 ], 0 );
					rr.move( 2*rayVector[ 1 ], 1 );
					rr.move( 2*rayVector[ 2 ], 2 );

					// refractive index at the projected location of the ray
					final double n1 = rr.get().get();

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

		return inject;
	}

	public static ImagePlus drawSingleRay( final RandomAccessibleInterval< FloatType > imgIn, final int z, final int x )
	{
		final ImageStack stack = new ImageStack( (int)imgIn.dimension( 0 ), (int)imgIn.dimension( 1 ) );

		// the refracted image
		final Img< FloatType > imgOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );
		final Img< FloatType > weightOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );

		final double[] sigma = new double[] { 0.5, 0.5, 0.5 };
		final VolumeInjection inject = new VolumeInjection( imgOut, weightOut, sigma );

		// the incoming direction of the light
		final double[] rayVector = new double[ 3 ];

		// the outgoing, refracted direction of the light
		final double[] refractedRay = new double[ 3 ];

		// the current position of the ray
		final double[] rayPosition = new double[ imgIn.numDimensions() ];

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		//final RandomAccess< FloatType > ra = img.randomAccess();
		final RealRandomAccess< FloatType > rr = Views.interpolate( Views.extendMirrorSingle( imgIn ), new NLinearInterpolatorFactory<>() ).realRandomAccess();

		final double valueIm = 1.0;

		rayVector[ 0 ] = 0;
		rayVector[ 1 ] = -1;
		rayVector[ 2 ] = 0;

		Raytrace.norm( rayVector );

		rayPosition[ 0 ] = x;
		rayPosition[ 1 ] = imgIn.dimension( 1 ) - 25;
		rayPosition[ 2 ] = z;

		while ( SimulateMultiViewAberrations.inside( rayPosition, imgIn ) )
		{
			rr.setPosition( rayPosition );

			// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)
			Hessian.computeHessianMatrix3D( rr, matrix );

			double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

			if ( Math.abs( ev ) > 0.01 )
			{
				// compute refractive index change
				rr.setPosition( rayPosition );

				rr.move( -rayVector[ 0 ], 0 );
				rr.move( -rayVector[ 1 ], 1 );
				rr.move( -rayVector[ 2 ], 2 );

				// refractive index at the origin of the ray
				final double n0 = rr.get().get();

				rr.move( 2*rayVector[ 0 ], 0 );
				rr.move( 2*rayVector[ 1 ], 1 );
				rr.move( 2*rayVector[ 2 ], 2 );

				// refractive index at the projected location of the ray
				final double n1 = rr.get().get();

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

			final FloatProcessor fp = getProcessor( inject.getWeight(), z );

			DecimalFormat df = new DecimalFormat("####0.00");
	
			Toolbar.setForegroundColor(java.awt.Color.WHITE);
			TextRoi tr = new TextRoi( 10,imgIn.dimension( 1 )-20,"vec=[" + df.format( eigenVector[ 0 ] ) + "," + df.format( eigenVector[ 0 ] ) + "], ev="+ df.format( ev ) );
			tr.setStrokeColor(Color.white); 
			tr.setNonScalable(true); 
			tr.drawPixels(fp);

			stack.addSlice( fp );

			System.out.println( Util.printCoordinates( rayPosition ) );
			rayPosition[ 0 ] += rayVector[ 0 ];
			rayPosition[ 1 ] += rayVector[ 1 ];
			rayPosition[ 2 ] += rayVector[ 2 ];
		}

		return new ImagePlus( "", stack );
	}

	public static VolumeInjection drawRays( final RandomAccessibleInterval< FloatType > imgIn, final int z, final int spacing )
	{
		// the refracted image
		final Img< FloatType > imgOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );
		final Img< FloatType > weightOut = new ArrayImgFactory< FloatType >( new FloatType() ).create( imgIn );

		final double[] sigma = new double[] { 0.5, 0.5, 0.5 };
		final VolumeInjection inject = new VolumeInjection( imgOut, weightOut, sigma );

		// the incoming direction of the light
		final double[] rayVector = new double[ 3 ];

		// the outgoing, refracted direction of the light
		final double[] refractedRay = new double[ 3 ];

		// the current position of the ray
		final double[] rayPosition = new double[ imgIn.numDimensions() ];

		final double[][] matrix = new double[ 3 ][ 3 ]; // row, column
		final double[] eigenVector = new double[ 3 ];

		//final RandomAccess< FloatType > ra = img.randomAccess();
		final RealRandomAccess< FloatType > rr = Views.interpolate( Views.extendMirrorSingle( imgIn ), new NLinearInterpolatorFactory<>() ).realRandomAccess();

		final double valueIm = 1.0;

		for ( int x = spacing; x <= imgIn.dimension( 0 ); x += spacing )
		{
			rayVector[ 0 ] = 0;
			rayVector[ 1 ] = -1;
			rayVector[ 2 ] = 0;

			Raytrace.norm( rayVector );

			rayPosition[ 0 ] = x;
			rayPosition[ 1 ] = imgIn.dimension( 1 ) - 25;
			rayPosition[ 2 ] = z;

			while ( SimulateMultiViewAberrations.inside( rayPosition, imgIn ) )
			{
				rr.setPosition( rayPosition );

				// normal vector of refraction plane (still maybe needs to be inverted to point towards the incoming signal)
				Hessian.computeHessianMatrix3D( rr, matrix );

				double ev = Hessian.computeLargestEigenVectorAndValue3d( matrix, eigenVector );

				if ( Math.abs( ev ) > 0.01 )
				{
					// compute refractive index change
					rr.setPosition( rayPosition );

					rr.move( -rayVector[ 0 ], 0 );
					rr.move( -rayVector[ 1 ], 1 );
					rr.move( -rayVector[ 2 ], 2 );

					// refractive index at the origin of the ray
					final double n0 = rr.get().get();

					rr.move( 2*rayVector[ 0 ], 0 );
					rr.move( 2*rayVector[ 1 ], 1 );
					rr.move( 2*rayVector[ 2 ], 2 );

					// refractive index at the projected location of the ray
					final double n1 = rr.get().get();

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

		return inject;
	}

	public static FloatProcessor getProcessor( final RandomAccessibleInterval< FloatType > img, final int z )
	{
		final RandomAccess< FloatType > ra = Views.zeroMin( img ).randomAccess();
		ra.setPosition( z, 2 );

		final FloatProcessor fp = new FloatProcessor( (int)img.dimension( 0 ), (int)img.dimension( 1 ));

		for ( int y = 0; y < img.dimension( 1 ); ++y )
		{
			ra.setPosition( y, 1 );

			for ( int x = 0; x < img.dimension( 0 ); ++x )
			{
				ra.setPosition( x, 0 );

				fp.setf( x, y, ra.get().get() );
			}
		}

		return fp;
	}

	public static void renderDifferentNAs()
	{
		final ImageStack stack = new ImageStack( 300, 300 );

		for ( double na = 1.0; na <= 2; na += 0.01 )
		{
			final Img< FloatType > img1 = ArrayImgs.floats( 300, 300, 10 );
			Raytrace.drawSimpleImage( img1, 1.0f, (float)na );

			final VolumeInjection rays = drawRays( img1, 5, 8 );

			final FloatProcessor fp = getProcessor( rays.getWeight(), 5 );

			Toolbar.setForegroundColor(java.awt.Color.WHITE);
			TextRoi tr = new TextRoi( 10,280,"NA="+ na );
			tr.setStrokeColor(Color.white); 
			tr.setNonScalable(true); 
			tr.drawPixels(fp);
			
			stack.addSlice( fp );
		}

		
		new ImagePlus( "", stack ).show();
	}

	public static void main( String[] args )
	{
		new ImageJ();

		//renderDifferentNAs();

		final Img< FloatType > img1 = ImageJFunctions.wrap( new ImagePlus( new File("spheres2b.tif").getAbsolutePath() ) );
		//final Img< FloatType > img1 = ArrayImgs.floats( 300, 300, 10 );
		//Raytrace.drawSimpleImage( img1, 1.0f, 1.33f );
		//drawSingleRay( img1, 150, 150 ).show();
		final VolumeInjection rays = drawMultipleRays( img1, 150, 150, 100 );
		ImageJFunctions.show( rays.getWeight() );

/*

		final Img< FloatType > img1 = ArrayImgs.floats( 300, 300, 300 );
		Raytrace.drawSimpleImage( img1, 1.0f, 1.33f );
		//ImageJFunctions.show( img1 );

		//final Img< FloatType > img2 = hessian( ImageJFunctions.wrap( new ImagePlus( new File("neuron.tif").getAbsolutePath() ) ), 0 );
		//final Img< FloatType > img2 = hessian( img1, 0 );
		//ImageJFunctions.show( img2 );

		final VolumeInjection rays = drawRays( img1, 5, 8 );
		ImageJFunctions.show( rays.getWeight() );*/
	}

}
