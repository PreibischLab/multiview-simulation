package simulation;

import java.util.ArrayList;
import java.util.Random;

import mpicbg.models.AffineModel3D;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

/**
 * Simulate Bead Images
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
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Stephan Preibisch (stephan.preibisch@gmx.de)
 */
public class SimulateBeads
{
	ArrayList< Img< FloatType > > imgs;
	public final int[] angles;
	public final int axis;
	public final int numPoints;
	public final Interval rangeSimulation;
	public final Interval intervalRender;
	public final double[] sigma;

	final Random rnd = new Random( 535 );

	public SimulateBeads( final int[] angles, final int axis, final int numPoints, final Interval rangeSimulation, final Interval intervalRender, final double[] sigma )
	{
		this.angles = angles;
		this.axis = axis;
		this.numPoints = numPoints;
		this.rangeSimulation = rangeSimulation;
		this.intervalRender = intervalRender;
		this.sigma = sigma;

		this.imgs = null;
	}

	public ArrayList< Img< FloatType > > getImgs()
	{
		if ( imgs == null )
		{
			final ArrayList< float[] > points = randomPoints( numPoints, rangeSimulation, rnd );
	
			final ArrayList< ArrayList< float[] > > lists = transformPoints( points, angles, axis, rangeSimulation );
	
			imgs = renderPoints( lists, intervalRender, sigma );
		}

		return imgs;
	}

	public static ArrayList< Img< FloatType > > renderPoints( final ArrayList< ArrayList< float[] > > lists, final Interval interval, final double[] sigma )
	{
		final ArrayList< Img< FloatType > > images = new ArrayList< Img< FloatType > >();

		for ( final ArrayList< float[] > list : lists )
		{
			final long[] dim = new long[ interval.numDimensions() ];

			for ( int d = 0; d < interval.numDimensions(); ++d )
				dim[ d ] = interval.max( d ) - interval.min( d );

			final Img< FloatType > img = ArrayImgs.floats( dim );

			for ( final float[] p : list )
				if ( isInsideAdjust( p, interval ) )
					addGaussian( img, p, sigma );

			images.add( img );
		}

		return images;
	}

	protected static boolean isInsideAdjust( final float[] p, final Interval interval )
	{
		for ( int d = 0; d < interval.numDimensions(); ++d )
		{
			p[ d ] -= interval.min( d );

			if ( p[ d ] < 0 || p[ d ] > interval.dimension( d ) - 1 )
				return false;
		}

		return true;
	}
	public static ArrayList< ArrayList< float[] > > transformPoints( final ArrayList< float[] > points, final int[] angles, final int axis, final Interval range )
	{
		final ArrayList< ArrayList< float[] > > transformedPoints = new ArrayList< ArrayList< float[] > >();

		for ( final int angle : angles )
		{
			final AffineModel3D t = SimulateMultiViewDataset.axisRotation( range, axis, angle );

			final ArrayList< float[] > transformed = new ArrayList< float[] >();

			for ( final float[] p : points )
				transformed.add( t.apply( p ) );

			transformedPoints.add( transformed );
		}

		return transformedPoints;
	}

	public static ArrayList< float[] > randomPoints( final int numPoints, final Interval range, final Random rnd )
	{
		final ArrayList< float[] > points = new ArrayList< float[] >();

		for ( int i = 0; i < numPoints; ++i )
		{
			final float[] p = new float[ range.numDimensions() ];

			for ( int d = 0; d < range.numDimensions(); ++d )
				p[ d ] = rnd.nextFloat() * ( range.max( d ) - range.min( d ) ) + range.min( d );

			points.add( p );
		}

		return points;
	}

	final public static void addGaussian( final Img< FloatType > image, final float[] location, final double[] sigma )
	{
		final int numDimensions = image.numDimensions();
		final int[] size = new int[ numDimensions ];
		
		final long[] min = new long[ numDimensions ];
		final long[] max = new long[ numDimensions ];
		
		final double[] two_sq_sigma = new double[ numDimensions ];
		
		for ( int d = 0; d < numDimensions; ++d )
		{
			size[ d ] = Util.getSuggestedKernelDiameter( sigma[ d ] ) * 2;
			min[ d ] = (int)Math.round( location[ d ] ) - size[ d ]/2;
			max[ d ] = min[ d ] + size[ d ] - 1;
			two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];
		}

		final RandomAccessible< FloatType > infinite = Views.extendZero( image );
		final RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
		final IterableInterval< FloatType > iterable = Views.iterable( interval );
		final Cursor< FloatType > cursor = iterable.localizingCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			
			double value = 1;
			
			for ( int d = 0; d < numDimensions; ++d )
			{
				final double x = location[ d ] - cursor.getIntPosition( d );
				value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
			}
			
			cursor.get().set( cursor.get().get() + ( (float)value * 1000.0f ) );
		}
	}

	public static void main( String[] args )
	{
		final int[] angles = new int[]{ 0, 45, 90, 135 };
		final int axis = 0;
		final int numPoints = 1000;
		final double[] sigma = new double[]{ 1, 1, 3 };
		final FinalInterval range = new FinalInterval( 512, 512, 200 );

		final SimulateBeads sb = new SimulateBeads( angles, axis, numPoints, range, range, sigma );

		for ( int i = 0; i < angles.length; ++i )
		{
			//ImageJFunctions.wrap( sb.getImgs().get( i ), "Angle_" + angles[ i ] ).duplicate().show();
			Tools.save( sb.getImgs().get( i ), "Angle_" + angles[ i ] + ".tif" );
		}

		System.out.println( "done." );
	}
}
