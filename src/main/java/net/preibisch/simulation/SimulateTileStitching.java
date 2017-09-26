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

import java.util.Random;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Pair;
import net.imglib2.util.Util;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

public class SimulateTileStitching
{
	final String dir = "src/main/resources/";

	int lightsheetSpacing = 3;
	float attenuation = 0.01f;

	final Random rnd;
	Img< FloatType > psf;
	Img< FloatType > con, conHalfPixel;
	long[] min, max;
	boolean halfPixelOffset;
	int[] overlap;

	public SimulateTileStitching( final boolean halfPixelOffset, final double[] overlapRatio )
	{
		this( null, halfPixelOffset, overlapRatio );
	}

	public SimulateTileStitching( final Random rnd, final boolean halfPixelOffset, final double[] overlapRatio )
	{
		if ( rnd == null )
			this.rnd = new Random( 464232194 );
		else
			this.rnd = rnd;

		this.psf = Tools.open( dir + "Angle0.tif", true );

		init( overlapRatio, halfPixelOffset );
	}

	public void init( final double[] overlapRatio, boolean halfPixelOffset )
	{
		this.halfPixelOffset = halfPixelOffset;

		// artificially rendered object based on which everything is computed,
		// including the ground-truth image, which is rotated once by the angle offset
		// so that it is realistic
		final int seed = rnd.nextInt(); // same result for half pixel and not

		final Thread[] threads = new Thread[ 2 ];

		threads[ 0 ] = new Thread( new Runnable()
		{
			public void run()
			{
				System.out.println( "Thread 1: simulate " );
				final Img< FloatType > groundTruth = SimulateMultiViewDataset.simulate( false, new Random( seed ) );
				final Img< FloatType > att = SimulateMultiViewDataset.attenuate3d( groundTruth, attenuation );
				con = SimulateMultiViewDataset.convolve( att, psf );
				Tools.adjustImage( con, SimulateMultiViewDataset.minValue, SimulateMultiViewDataset.avgIntensity );
			}
		} );

		threads[ 1 ] = new Thread( new Runnable()
		{
			public void run()
			{
				System.out.println( "Thread 2: simulate " );
				final Img< FloatType > groundTruthHalfPixel = SimulateMultiViewDataset.simulate( true, new Random( seed ) );
				final Img< FloatType > attHalfPixel = SimulateMultiViewDataset.attenuate3d( groundTruthHalfPixel, attenuation );
				conHalfPixel = SimulateMultiViewDataset.convolve( attHalfPixel, psf );
				Tools.adjustImage( conHalfPixel, SimulateMultiViewDataset.minValue, SimulateMultiViewDataset.avgIntensity );
			}
		} );

		runThreads( threads );

		this.overlap = new int[ con.numDimensions() ];
		for ( int d = 0; d < con.numDimensions(); ++d )
			this.overlap[ d ] = (int)Math.round( con.dimension( d ) * overlapRatio[ d ] / 2 );

		System.out.println( "Overlap: " + Util.printCoordinates( this.overlap ) );

		this.min = new long[ con.numDimensions() ];
		this.max = new long[ con.numDimensions() ];
	}

	private Img< FloatType > split0, split1;

	public Pair< Img< FloatType >, Img< FloatType > > getNextPair( final float snr )
	{
		final Thread[] threads = new Thread[ 2 ];

		final int seed0 = rnd.nextInt();
		final int seed1 = rnd.nextInt();

		// tile 0
		threads[ 0 ] = new Thread( new Runnable()
		{
			public void run()
			{
				final long[] min1 = min.clone();
				final long[] max1 = max.clone();

				getInterval( min1, max1, 0 );

				//System.out.println( new Date( System.currentTimeMillis() ) + ": extracting slices tile " + tile + " snr=" +snr );
				split0 = SimulateMultiViewDataset.extractSlices(
						Views.zeroMin( Views.interval( con, min1, max1 ) ),
						lightsheetSpacing,
						snr,
						new Random( seed0 ) );
			}
		} );

		// tile 1
		threads[ 1 ] = new Thread( new Runnable()
		{
			public void run()
			{
				final long[] min1 = min.clone();
				final long[] max1 = max.clone();

				getInterval( min1, max1, 1 );

				if ( halfPixelOffset )
					split1 = SimulateMultiViewDataset.extractSlices(
							Views.zeroMin( Views.interval( conHalfPixel, min1, max1 ) ),
							lightsheetSpacing, snr, new Random( seed1 ) );
				else
					split1 = SimulateMultiViewDataset.extractSlices(
							Views.zeroMin( Views.interval( con, min1, max1 ) ),
							lightsheetSpacing, snr, new Random( seed1 ) );
			}
		} );

		runThreads( threads );

		return new ValuePair<>( split0, split1 );
	}

	public double[] getCorrectTranslation()
	{
		final double[] translation = new double[ con.numDimensions() ];

		getInterval( min, max, 1 );

		for ( int d = 0; d < con.numDimensions(); ++d )
			translation[ d ] = min[ d ];

		if ( halfPixelOffset )
		{
			// unintutive, but true. If the content of the right image is shifted
			// more to the right (by 0.5 pixels), means that in order to match it
			// has to be shifted less to the left (take to sheets of paper to understand :)
			translation[ 0 ] -= 0.5;
			translation[ 1 ] -= 0.5;
		}

		translation[ 2 ] /= lightsheetSpacing;

		return translation;
	}

	protected void getInterval( final long[] min, final long[] max, final int tile )
	{
		con.min( min );
		con.max( max );

		if ( tile == 0 )
		{
			for ( int d = 0; d < min.length; ++d )
				max[ d ] = (int)con.dimension( d ) / 2 + overlap[ d ];

			System.out.println( "Tile0: " + Util.printCoordinates( min ) + " >> " + Util.printCoordinates( max ) );
		}
		else
		{
			for ( int d = 0; d < min.length; ++d )
				min[ d ] = (int)con.dimension( 0 ) / 2 - overlap[ d ];

			System.out.println( "Tile1: " + Util.printCoordinates( min ) + " >> " + Util.printCoordinates( max ) );
		}

	}

	public static void show( final RandomAccessibleInterval< FloatType > img, final String title )
	{
		ImagePlus imp = ImageJFunctions.wrapFloat( img, title );

		imp.setDimensions( 1, imp.getStackSize(), 1 );
		imp.setSlice( imp.getStackSize() / 2 );
		imp.resetDisplayRange();
		imp.show();
	}

	public static void runThreads( final Thread[] threads )
	{
		for ( int ithread = 0; ithread < threads.length; ++ithread )
			threads[ ithread ].start();

		try
		{
			for ( int ithread = 0; ithread < threads.length; ++ithread )
				threads[ ithread ].join();
		}
		catch ( InterruptedException ie ) { throw new RuntimeException(ie); }
	}

	public static void main( String[] args )
	{
		new ImageJ();

		final double overlap = 0.2;
		final float snr = 8;

		final SimulateTileStitching sts = new SimulateTileStitching( new Random( System.currentTimeMillis() ), true, Util.getArrayFromValue( overlap, 3 ) );
		
		IJ.log( "Known shift (right relative to left): " + Util.printCoordinates( sts.getCorrectTranslation() ) );

		final Pair< Img< FloatType >, Img< FloatType > > pair = sts.getNextPair( snr );

		show( pair.getA(), "left" );
		show( pair.getB(), "right" );
	}
}
