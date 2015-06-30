package simulation;

import ij.ImagePlus;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.awt.Image;
import java.util.Random;

import org.uncommons.maths.random.PoissonGenerator;

import mpicbg.util.RealSum;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

/**
 * static helper methods
 * @author Stephan Preibisch (stephan.preibisch@gmx.de)
 *
 */
public class Tools
{
	public static void poissonProcess( final RandomAccessibleInterval<FloatType> img, final double SNR, final Random rnd )
	{
		// based on an average intensity of 5, a multiplicator of 1 corresponds to a SNR of 2.23 = sqrt( 5 );	
		final double mul = Math.pow( SNR / Math.sqrt( 5 ), 2 );
		
		final NumberGeneratorImage< FloatType> ng = new NumberGeneratorImage< FloatType>( img, mul );
		final PoissonGenerator pg = new PoissonGenerator( ng, rnd );
		
		for ( final FloatType v : Views.iterable( img ) )
		{
			ng.fwd();
			v.set( pg.nextValue().floatValue() );
		}
	}

    public static void save( final Img< FloatType > img, final String file )
    {
    	final ImagePlus imp = ImageJFunctions.wrap( img, "" ).duplicate();

    	final FileSaver s = new FileSaver( imp );
    	
    	if ( img.numDimensions() == 3 )
    	{
    		imp.setDimensions( 1, (int)img.dimension( 2 ), 1 );
    		s.saveAsTiffStack( file );
    	}
    	else
    	{
    		s.saveAsTiff( file );
    	}

    	imp.close();
    }
	
	/**
	 * Norms an image so that the sum over all pixels is 1.
	 * 
	 * @param img - the {@link Image} to normalize
	 */
	final public static void normImage( final Iterable<FloatType> img )
	{
		final double sum = sumImage( img );	

		for ( final FloatType t : img )
			t.set( (float) ((double)t.get() / sum) );
	}

	/**
	 * @param img - the input {@link Image}
	 * @return - the sum of all pixels using {@link RealSum}
	 */
	final public static double sumImage( final Iterable<FloatType> img )
	{
		final RealSum sum = new RealSum();		

		for ( final FloatType t : img )
			sum.add( t.get() );

		return sum.getSum();
	}	

	/**
	 * Adjusts an image so that the minimal intensity is minValue and the average is average
	 * 
	 * @param image - the image to norm
	 * @param minValue - the minimal value
	 * @param targetAverage - the average that we want to have
	 * 
	 * @return by which number all intensities were multiplied
	 */
	public static double adjustImage( final IterableInterval<FloatType> image, final float minValue, final float targetAverage )
	{
		// first norm the image to an average of (targetAverage - minValue)
		final double avg = sumImage( image )/(double)image.size();
		final double correction = ( targetAverage - minValue ) / avg;

		// correct 
		for ( final FloatType t : image )
			t.set( (float)( t.get() * correction ) );
			
		// now add minValue to all pixels
		for ( final FloatType t : image )
			t.set( t.get() + minValue );
		
		//System.out.println( sumImage( image )/(double)image.getNumPixels() );
		return correction;
	}

	public static Img< FloatType > open( final String name, final ImgFactory< FloatType > factory )
	{
		final Opener io = new Opener();
		ImagePlus imp = io.openImage( name );
		
		if ( imp.getStack().getSize() > 1 )
		{
			final int depth = imp.getStack().getSize();
			
			final Img<FloatType> img = factory.create( new int[]{ imp.getWidth(), imp.getHeight(), depth }, new FloatType() );
			
			final RandomAccess< FloatType > r = img.randomAccess();
			final int[] tmp = new int[]{ 0, 0, 0 };
			
			for ( int i = 0; i < depth; ++i )
			{
				final ImageProcessor ip = imp.getStack().getProcessor( i + 1 );
				tmp[ 2 ] = i;
				tmp[ 0 ] = tmp[ 1 ] = 0;
				
				r.setPosition( tmp );
				
				for ( int y = 0; y < imp.getHeight(); ++y )
				{
					r.setPosition( y, 1 );
					for ( int x = 0; x < imp.getWidth(); ++x )
					{
						r.setPosition( x, 0 );
						r.get().set( ip.getPixelValue( x, y ) );
					}
				}
			}			

			return img;
		}
		else
		{
			final Img<FloatType> img = factory.create( new int[]{ imp.getWidth(), imp.getHeight() }, new FloatType() );
			
			final ImageProcessor ip = imp.getProcessor();
			
			final Cursor<FloatType> cursorOut = img.cursor();
			
			while ( cursorOut.hasNext() )
			{
				cursorOut.fwd();
				cursorOut.get().set( ip.getPixelValue( cursorOut.getIntPosition( 0 ), cursorOut.getIntPosition( 1 ) ) );
			}
			
			return img;
		}
	}

    public static Img< FloatType > open( final String file, final boolean square )
    {
    	if ( square )
    	{
    		return Tools.makeSquare( open( file, new ArrayImgFactory<FloatType>() ) );
    	}
    	else
    	{
    		return open( file, new ArrayImgFactory<FloatType>() );
    	}
    }

	/**
	 * Make quadratic image, pad with minimal intensity
	 * 
	 * @param img
	 * @return
	 */
	public static Img< FloatType > makeSquare( final Img< FloatType > img )
	{
		final long[] tmp = new long[ img.numDimensions() ];
		long maxSize = 0;
		
		for ( int d = 0; d < img.numDimensions(); ++d )
			maxSize = Math.max( maxSize, img.dimension( d ) );
		
		for ( int d = 0; d < img.numDimensions(); ++d )
			tmp[ d ] = maxSize;
		
		float min = Float.MAX_VALUE;
		
		for ( final FloatType f : img )
			min = Math.min( min, f.get() );
		
		final Img< FloatType > square = img.factory().create( tmp, img.firstElement() );
		
		final Cursor< FloatType > squareCursor = square.localizingCursor();
		final RandomAccess< FloatType > inputCursor = Views.extendValue( img, new FloatType( min ) ).randomAccess(); 
		
		while ( squareCursor.hasNext() )
		{
			squareCursor.fwd();
			squareCursor.localize( tmp );
			
			for ( int d = 0; d < img.numDimensions(); ++d )
				tmp[ d ] =  tmp[ d ] - square.dimension( d )/2 + img.dimension( d )/2;

			inputCursor.setPosition( tmp );
			squareCursor.get().set( inputCursor.get().get() );
		}
		
		return square;
	}

}
