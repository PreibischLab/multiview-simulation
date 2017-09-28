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
package net.preibisch.simulation.imgloader;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import mpicbg.spim.data.legacy.LegacyImgLoader;
import mpicbg.spim.data.sequence.FinalVoxelDimensions;
import mpicbg.spim.data.sequence.ViewId;
import mpicbg.spim.data.sequence.VoxelDimensions;
import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.type.numeric.real.FloatType;
import net.preibisch.simulation.SimulateBeads;
import net.preibisch.simulation.SimulateBeads2;

public class LegacySimulatedBeadsImgLoader2 implements LegacyImgLoader< UnsignedShortType >
{
	final SimulateBeads2 sb;
	final Map<ViewId, List<Integer>> vidToDescription;
	
	public LegacySimulatedBeadsImgLoader2( final SimulateBeads2 sb ) 
	{
		this.sb = sb;
		this.vidToDescription = new HashMap<>();
	}
	
	
	public void addViewId(int idx, int tp, int angle, int channel, int tile, int illumination)
	{
		List< Integer > val = Arrays.asList( new Integer[] {tp, angle, channel, tile, illumination} );
		vidToDescription.put( new ViewId(tp, idx),  val );
	}

	public SimulateBeads2 getSimulateBeads() { return sb; }

	@Override
	public RandomAccessibleInterval< UnsignedShortType > getImage(ViewId view)
	{
		List< Integer > key = vidToDescription.get( view );
		System.out.println( key );
		Img< FloatType > imgF = sb.getImg(key.get( 0 ), key.get( 1 ), key.get( 2 ), key.get( 3 ), key.get( 4 ));
		final long[] dim = new long[ imgF.numDimensions() ];

		for ( int d = 0; d < dim.length; ++d )
			dim[ d ] = imgF.dimension( d );

		final Img< UnsignedShortType > img = ArrayImgs.unsignedShorts( dim );

		final Cursor< FloatType > in = imgF.cursor();
		final Cursor< UnsignedShortType > out = img.cursor();

		while ( in.hasNext() )
			out.next().set( Math.round( in.next().get() ) );

		return img;
	}

	@Override
	public UnsignedShortType getImageType()
	{
		return new UnsignedShortType();
	}
	

	@Override
	public RandomAccessibleInterval< FloatType > getFloatImage(ViewId view, boolean normalize)
	{
		List< Integer > key = vidToDescription.get( view );
		Img< FloatType > imgF = sb.getImg(key.get( 0 ), key.get( 1 ), key.get( 2 ), key.get( 3 ), key.get( 4 ));
		if ( normalize )
		{
			return LegacySimulatedBeadsImgLoader.normalize( imgF.copy() );
		}
		else
		{
			return imgF.copy();
		}
	}

	@Override
	public Dimensions getImageSize(ViewId view)
	{
		List< Integer > key = vidToDescription.get( view );
		Img< FloatType > imgF = sb.getImg(key.get( 0 ), key.get( 1 ), key.get( 2 ), key.get( 3 ), key.get( 4 ));
		return imgF;
	}

	@Override
	public VoxelDimensions getVoxelSize(ViewId view)
	{
		return new FinalVoxelDimensions( "pixel", 1, 1, 1 );
	}
	
	
	

}
