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

import java.awt.Image;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.ImageProcessor;
import loci.formats.ChannelSeparator;
import loci.formats.FormatException;
import loci.formats.FormatTools;
import loci.formats.IFormatReader;
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
import net.preibisch.simulation.uncommons.PoissonGenerator;

/**
 * static helper methods
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

    public static void save( final RandomAccessibleInterval< FloatType > img, final String file )
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
		final IFormatReader r = new ChannelSeparator();

		final String id = new File( name ).getAbsolutePath();

		try
		{
			r.setId( id );
			
			final boolean isLittleEndian = r.isLittleEndian();
			final int width = r.getSizeX();
			final int height = r.getSizeY();
			final int depth = r.getSizeZ();
			int timepoints = r.getSizeT();
			int channels = r.getSizeC();
			final int tiles = r.getSeriesCount();
			final int pixelType = r.getPixelType();
			final int bytesPerPixel = FormatTools.getBytesPerPixel( pixelType );
			final String pixelTypeString = FormatTools.getPixelTypeString( pixelType );

			if ( pixelType != FormatTools.FLOAT )
			{
				System.out.println( "StackImgLoaderLOCI.openLOCI(): PixelType " + pixelTypeString + " not supported, returning. ");

				r.close();

				return null;
			}

			System.out.println( "w: " + width );
			System.out.println( "h: " + height );
			System.out.println( "d: " + depth );
			System.out.println( "t: " + timepoints );
			System.out.println( "c: " + channels );

			final Img< FloatType > img = factory.create( new long[] { width, height, depth }, new FloatType() );

			final byte[] b = new byte[width * height * bytesPerPixel];

			final int planeX = 0;
			final int planeY = 1;

			for ( int z = 0; z < depth; ++z )
			{
				final Cursor< FloatType > cursor = Views.iterable( Views.hyperSlice( img, 2, z ) ).localizingCursor();

				r.openBytes( r.getIndex( z, 0, 0 ), b );

				while( cursor.hasNext() )
				{
					cursor.fwd();
					cursor.get().setReal( getFloatValue( b, ( cursor.getIntPosition( planeX )+ cursor.getIntPosition( planeY )*width )*4, isLittleEndian ) );
				}
			}

			r.close();

			return img;

		} catch ( FormatException | IOException e )
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		return null;
	}

	public static final float getFloatValue( final byte[] b, final int i, final boolean isLittleEndian )
	{
		if ( isLittleEndian )
			return Float.intBitsToFloat( ((b[i+3] & 0xff) << 24)  + ((b[i+2] & 0xff) << 16)  +  ((b[i+1] & 0xff) << 8)  + (b[i] & 0xff) );
		else
			return Float.intBitsToFloat( ((b[i] & 0xff) << 24)  + ((b[i+1] & 0xff) << 16)  +  ((b[i+2] & 0xff) << 8)  + (b[i+3] & 0xff) );
	}


	public static Img< FloatType > openIJ( final String name, final ImgFactory< FloatType > factory )
	{
		if ( !new File( name ).exists() )
			throw new RuntimeException( "file '" + name + "' does not exisit." );

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
	 * @param img - the input
	 * @return - the squared image
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
