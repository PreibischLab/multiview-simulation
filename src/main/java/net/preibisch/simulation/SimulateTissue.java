package net.preibisch.simulation;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.RandomAccessibleOnRealRandomAccessible;
import net.imglib2.view.Views;

public class SimulateTissue
{
	private static class SimulationParameters
	{
		long seed = 42;
		long[] max = new long[]{1024, 1024, 256};
		long[] min = new long[]{0, 0, 0};
		
		double[] perlinScales = new double[] {1024/4,1025/1.5,256};
		double perlinThresh = 0.15;
		double backgroundDensity = 0.01;

		int nBigSpheres = 200;
		float minRadiusBig = 20;
		float maxRadiusBig = 35;
		float minValueBig = 1.2f;
		float maxValueBig = 2.4f;

		int nSmallSamples = 20000;
		float minRadiusSmall = 2;
		float maxRadiusSmall = 4;
		float minValueSmall = 4.0f;
		float maxValueSmall = 6.0f;

		double[] transform = new AffineTransform3D().getRowPackedCopy();

		public SimulationParameters(){}
	}

	public static void main(String[] args)
	{
		// load from JSON if given
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		SimulationParameters paramsT = null;
		if (args.length < 1) { return; }
		else
		{
			try (FileReader fr = new FileReader( new File(args[0], "sim_params.json").getPath() ))
			{
				paramsT = gson.fromJson( fr, SimulationParameters.class);
			}
			catch ( FileNotFoundException e )
			{
				e.printStackTrace();
			}
			catch ( IOException e )
			{
				e.printStackTrace();
			}
		}
		final SimulationParameters params = paramsT;

		// create out directory
		if (!Files.exists( Paths.get( args[0], "sim-phantom" )))
			try
			{
				Files.createDirectory( Paths.get( args[0], "sim-phantom" ));
			}
			catch ( IOException e )
			{
				e.printStackTrace();
			}

		// seeded RNG for reproducuibility
		final Random rnd = new Random( params.seed );

		// general tissue shape: thresholded Perlin noise
		RealRandomAccessible< FloatType > rrablePerlin = new PerlinNoiseRealRandomAccessible<>( new FloatType(), params.perlinScales, new int[] {15, 15, 15}, 100, rnd );
		RealRandomAccessible< FloatType > rrablePerlinThrd = new SimpleCalculatedRealRandomAccessible<FloatType>( new FloatType(), (a,b) -> {
			a.setReal( b.iterator().next().get() > params.perlinThresh ? 1.0 : 0);
		}, rrablePerlin);

		// create big spherical objects within tissue
		// these serve as density for sampling of smaller points
		List< RealPoint > bigSpherePositions = PointRejectionSampling.sampleRealPoints( new FinalInterval( params.min, params.max ), params.nBigSpheres, rrablePerlinThrd, rnd );
		List<Double> bigSphereRadii = new ArrayList<>();
		HypersphereCollectionRealRandomAccessible< FloatType > rrableDensity = new HypersphereCollectionRealRandomAccessible<>( params.min.length, new FloatType() );
		for (int i=0; i<params.nBigSpheres; i++)
		{
			double radius = params.minRadiusBig + rnd.nextDouble() * (params.maxRadiusBig - params.minRadiusBig);
			rrableDensity.addSphere( 
					bigSpherePositions.get( i ),
					radius,
					new FloatType(1.0f) );
			bigSphereRadii.add( radius );
		}

		// add some density to background -> smaller objects will no be exclusively within bigger objs
		RealRandomAccessible< FloatType > rrableDensity2 = new SimpleCalculatedRealRandomAccessible<FloatType>( new FloatType(), (a,b) -> {
			a.setReal( b.iterator().next().get() > 0.0 ? 1.0 : params.backgroundDensity);
		}, rrableDensity);

		HypersphereCollectionRealRandomAccessible< FloatType > rrableBigPoints = new HypersphereCollectionRealRandomAccessible<>( params.min.length, new FloatType() );
		HypersphereCollectionRealRandomAccessible< FloatType > rrableSmallPoints = new HypersphereCollectionRealRandomAccessible<>( params.min.length, new FloatType() );
		List< RealPoint > smallPoints = PointRejectionSampling.sampleRealPoints( new FinalInterval( params.min, params.max ), params.nSmallSamples, rrableDensity2, rnd );
		// create separate rrables for big and small objects
		// (separate for performance reasons: HypersphereCollectionRealRandomAccessible is fast for many small objects,
		//  but not for many small and some big objects)
		for (final RealPoint sp : smallPoints)
		{
			rrableSmallPoints.addSphere( 
					sp,
					params.minRadiusSmall + rnd.nextDouble() * (params.maxRadiusSmall - params.minRadiusSmall),
					new FloatType((float) ( params.minValueSmall + rnd.nextDouble() * (params.maxValueSmall - params.minValueSmall) )) );
		}
		for (int i=0; i<params.nBigSpheres; i++)
		{
			rrableBigPoints.addSphere( 
					bigSpherePositions.get( i ),
					bigSphereRadii.get( i ),
					new FloatType((float) ( params.minValueBig + rnd.nextDouble() * (params.maxValueBig - params.minValueBig) )) );
		}

		// final image: maximum of tissue img, big objects, small objects imgs
		RealRandomAccessible< FloatType > rrableFinalUntransformed = new SimpleCalculatedRealRandomAccessible<FloatType>( new FloatType(), (a,b) -> {
			float res = 0;
			for (FloatType t: b )
				res = Math.max( res, t.getRealFloat() );
			a.setReal( res );
		}, rrableBigPoints, rrableSmallPoints, rrablePerlinThrd);

		// transform the image with an arbitrary affine transform
		AffineTransform transform = new AffineTransform( params.min.length );
		transform.set( params.transform );
		final RealRandomAccessible< FloatType > rrableFinal = RealViews.affineReal( rrableFinalUntransformed, transform );

		RandomAccessibleOnRealRandomAccessible< FloatType > raster = Views.raster( rrableFinal );
		if (args.length > 1 && args[1].equals( "display" ))
		{
			// raster interval and display
			new ImageJ();
			ImageJFunctions.show( Views.interval( raster, new FinalInterval( params.min, params.max )),
					Executors.newFixedThreadPool( Runtime.getRuntime().availableProcessors() ) );
		}
		

		// save result as multiple tiffs (downstream python has problems reading 4gb+ tiff stacks)
		if (args.length > 1 && args[1].equals( "save" ))
		{
			ExecutorService service = Executors.newFixedThreadPool( Runtime.getRuntime().availableProcessors() );
			ArrayList<Callable< Void >> calls = new ArrayList<>();
			for (long z=0; z<(params.max[2] - params.min[2]); z++)
			{
				final long theZ = z;
				calls.add( new Callable< Void >()
				{
					@Override
					public Void call() throws Exception
					{
						// this may be a little reinvention of the wheel, but weird errors happened with Views.hyperSlice
						// raster interval and display
						final RealRandomAccess< FloatType > rra = rrableFinal.realRandomAccess();

						long[] minI = params.min.clone();
						long[] maxI = params.max.clone();

						minI[2] += theZ;
						maxI[2] = minI[2];

						RandomAccessibleInterval< FloatType> out = new ArrayImgFactory<>( new FloatType() ).create( new FinalInterval( minI, maxI ) );
						out = Views.translate( out, minI );
						
						Cursor< FloatType > c = Views.iterable( out ).localizingCursor();
						while(c.hasNext()) {
							c.fwd();
							rra.setPosition( c );
							c.get().set( rra.get() );
						}

						RandomAccessibleInterval< FloatType > slice = Views.dropSingletonDimensions( out );
						ImagePlus imp = ImageJFunctions.wrap(slice, "plane" + theZ);
						IJ.saveAsTiff( imp, Paths.get( args[0], "sim-phantom", "sim"+theZ+".tif" ).toString() );

						System.out.println( "saved plane " + (theZ+1) + " (of " + (params.max[2] - params.min[2]) + " planes)." );
						return null;
					}
				} );
			}
			List< Future< Void > > fs = null;
			try
			{
				fs = service.invokeAll( calls );
				for (Future< Void > f : fs)
					f.get();
			}
			catch ( InterruptedException | ExecutionException e )
			{
				e.printStackTrace();
			}
			System.out.println( "DONE." );
			System.exit( 0 );
		}


	}
}
