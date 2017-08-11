package simulation;

import java.util.ArrayList;
import java.util.Date;
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
	final Img< FloatType > con, conHalfPixel;
	final long[] min, max;
	final boolean halfPixelOffset;
	final int overlap;

	public SimulateTileStitching( final boolean halfPixelOffset )
	{
		this( null, halfPixelOffset, 0.2 );
	}

	public SimulateTileStitching( final boolean halfPixelOffset, final double overlapRatio )
	{
		this( null, halfPixelOffset, overlapRatio );
	}

	public SimulateTileStitching( final Random rnd, final boolean halfPixelOffset )
	{
		this( rnd, halfPixelOffset, 0.2 );
	}

	public SimulateTileStitching( final Random rnd, final boolean halfPixelOffset, final double overlapRatio )
	{
		this.halfPixelOffset = halfPixelOffset;

		if ( rnd == null )
			this.rnd = new Random( 464232194 );
		else
			this.rnd = rnd;

		// artificially rendered object based on which everything is computed,
		// including the ground-truth image, which is rotated once by the angle offset
		// so that it is realistic
		System.out.println( new Date( System.currentTimeMillis() ) + ": rendering basis for ground truth" );
		final Img< FloatType > groundTruth = SimulateMultiViewDataset.simulate( false, new Random( 464232194 ) );
		final Img< FloatType > groundTruthHalfPixel = SimulateMultiViewDataset.simulate( true, new Random( 464232194 ) );

		this.overlap = (int)Math.round( groundTruth.dimension( 0 ) * overlapRatio );
		System.out.println( "Overlap: " + overlap );

		this.min = new long[ groundTruth.numDimensions() ];
		this.max = new long[ groundTruth.numDimensions() ];

		System.out.println( new Date( System.currentTimeMillis() ) + ": attenuation " );
		Img<FloatType> att = SimulateMultiViewDataset.attenuate3d( groundTruth, attenuation );
		Img<FloatType> attHalfPixel = SimulateMultiViewDataset.attenuate3d( groundTruthHalfPixel, attenuation );

		//ImageJFunctions.show( att ).setTitle( "att" );

		System.out.println( new Date( System.currentTimeMillis() ) + ": convolving "  );
		Img< FloatType > psf = Tools.open( dir + "Angle0.tif", true );
		this.con = SimulateMultiViewDataset.convolve( att, psf );
		this.conHalfPixel = SimulateMultiViewDataset.convolve( attHalfPixel, psf );

		Tools.adjustImage( con, SimulateMultiViewDataset.minValue, SimulateMultiViewDataset.avgIntensity );
		Tools.adjustImage( conHalfPixel, SimulateMultiViewDataset.minValue, SimulateMultiViewDataset.avgIntensity );
	}

	public Pair< Img< FloatType >, Img< FloatType > > getNextPair( final float snr )
	{
		final ArrayList< Img< FloatType > > imgs = new ArrayList<>();

		for ( int tile = 0; tile <= 1; ++tile )
		{
			getInterval( min, max, tile );

			final RandomAccessibleInterval< FloatType > split;

			if ( tile == 1 && halfPixelOffset )
			{
				//System.out.println( new Date( System.currentTimeMillis() ) + ": extracting slices tile " + tile + " (0.5px offset), snr=" +snr );
				split = Views.zeroMin( Views.interval( conHalfPixel, min, max ) );
			}
			else
			{
				//System.out.println( new Date( System.currentTimeMillis() ) + ": extracting slices tile " + tile + " snr=" +snr );
				split = Views.zeroMin( Views.interval( con, min, max ) );
			}

			imgs.add( SimulateMultiViewDataset.extractSlices( split, lightsheetSpacing, snr ) );
		}

		return new ValuePair<>( imgs.get( 0 ), imgs.get( 1 ) );
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

		return translation;
	}

	protected void getInterval( final long[] min, final long[] max, final int tile )
	{
		con.min( min );
		con.max( max );

		if ( tile == 0 )
			max[ 0 ] = con.dimension( 0 ) / 2 + overlap;
		else
			min[ 0 ] = con.dimension( 0 ) / 2 - overlap;
	}

	public static void show( final RandomAccessibleInterval< FloatType > img, final String title )
	{
		ImagePlus imp = ImageJFunctions.wrapFloat( img, title );

		imp.setDimensions( 1, imp.getStackSize(), 1 );
		imp.setSlice( imp.getStackSize() / 2 );
		imp.resetDisplayRange();
		imp.show();
	}

	public static void main( String[] args )
	{
		new ImageJ();

		final double overlap = 0.2;
		final float snr = 8;

		final SimulateTileStitching sts = new SimulateTileStitching( new Random( System.currentTimeMillis() ), true, overlap );
		
		IJ.log( "Known shift (right relative to left): " + Util.printCoordinates( sts.getCorrectTranslation() ) );

		final Pair< Img< FloatType >, Img< FloatType > > pair = sts.getNextPair( snr );

		show( pair.getA(), "left" );
		show( pair.getB(), "right" );
	}
}
