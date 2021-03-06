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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.imglib2.Dimensions;
import net.imglib2.FinalRealInterval;
import net.imglib2.Interval;
import net.imglib2.RealInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;

public class SimulateBeads2
{

	public final Map<List<Integer>, Img< FloatType > > imgs;
	public final int numPoints;
	public final Interval rangeSimulation;
	public final Interval intervalRender;
	public final double[] sigma;
	
	public final Map<Integer, AffineTransform3D> angleTransforms;
	public final Map<Integer, AffineTransform3D> channelTransforms;
	public final Map<Integer, AffineTransform3D> illumTransforms;
	public final Map<Integer, AffineTransform3D> tpTransforms;
	public final Map<Integer, AffineTransform3D> tileTransforms;

	final Random rnd = new Random( 535 );
	final ArrayList< double[] > points;
	
	public SimulateBeads2(int numPoints, double[] sigma, Interval rangeSimulation, Interval intervalRender)
	{
		this.numPoints = numPoints;
		this.rangeSimulation = rangeSimulation;
		this.intervalRender = intervalRender;
		this.sigma = sigma;
				
		imgs = new HashMap<>();
		
		angleTransforms = new HashMap<>();
		channelTransforms = new HashMap<>();
		illumTransforms = new HashMap<>();
		tpTransforms = new HashMap<>();
		tileTransforms = new HashMap<>();
		
		points = SimulateBeads.randomPoints( numPoints, rangeSimulation, rnd );
	}
	
	// get min/max of all tiles
	public RealInterval getTilesExtent()
	{
		final Dimensions dims = intervalRender;
		double[] mins = Util.getArrayFromValue( Double.MAX_VALUE, dims.numDimensions() );
		double[] maxs = Util.getArrayFromValue( -Double.MAX_VALUE, dims.numDimensions() );
		
		for (AffineTransform3D tt : tileTransforms.values())
		{
			for (int d = 0; d < dims.numDimensions(); d++)
			{
				mins[d] = Math.min( mins[d], -tt.getTranslation()[d] );
				maxs[d] = Math.max( maxs[d], -tt.getTranslation()[d] + dims.dimension( d ) );
			}
		}
		
		return new FinalRealInterval( mins, maxs );
		
	}
	
	public Img< FloatType >  getImg(int tp, int angle, int channel, int tile, int illumination)
	{
		List< Integer > key = Arrays.asList( new Integer[] {tp, angle, channel, tile, illumination} );
		if (!imgs.containsKey( key ))
		{
			ArrayList< double[] > transformedPoints = transformedPoints( tp, angle, channel, tile, illumination );
			ArrayList< ArrayList< double[] > > packed = new ArrayList<>();
			packed.add( transformedPoints );
			imgs.put( key, SimulateBeads.renderPoints( packed, intervalRender, sigma ).get( 0 ));
		}
		return imgs.get( key );
		
	}
	
	public void addAngle(int id, int axis, double degrees)
	{
		AffineTransform3D t = new AffineTransform3D();
		t.rotate( axis, Math.toRadians( degrees ) );
		angleTransforms.put( id, t.copy() );
	}
	
	public void addChannel(int id, double[] shift)
	{
		AffineTransform3D t = new AffineTransform3D();
		t.translate( shift );
		channelTransforms.put( id, t.copy());
	}
	
	public void addIllumination(int id, double[] shift)
	{
		AffineTransform3D t = new AffineTransform3D();
		t.translate( shift );
		illumTransforms.put( id, t.copy() );
	}
	
	public void addTimepoint(int id, double[] shift)
	{
		AffineTransform3D t = new AffineTransform3D();
		t.translate( shift );
		tpTransforms.put( id, t.copy());
	}
	
	public void addTile(int id, double[] shift)
	{
		AffineTransform3D t = new AffineTransform3D();
		t.translate( shift );
		tileTransforms.put( id, t.copy().inverse() );
	}
	
	private ArrayList< double[] > transformedPoints(int tp, int angle, int channel, int tile, int illumination)
	{
		
		AffineTransform3D transform = new AffineTransform3D();
		
		// time
		if (tpTransforms.containsKey( tp ))
			transform.preConcatenate( tpTransforms.get( tp ) );
		
		// angle
		if (angleTransforms.containsKey( angle ))
			transform.preConcatenate( angleTransforms.get( angle ) );
		
		// channel
		if (channelTransforms.containsKey( channel ))
			transform.preConcatenate( channelTransforms.get( channel ) );
		
		// illumination
		if (illumTransforms.containsKey( illumination ))
			transform.preConcatenate( illumTransforms.get( illumination ) );
		
		// tile
		if (tileTransforms.containsKey( tile ))
			transform.preConcatenate( tileTransforms.get( tile ) );
				
		final ArrayList< double[] > transformed = new ArrayList< double[] >();

		for ( final double[] p : points )
		{
			double[] pt = p.clone();
			transform.apply( p, pt );
			transformed.add( pt );
		}
			
		return transformed;
		
	}
	
	public static void main(String[] args)
	{
		SimulateBeads2 sb = new SimulateBeads2( 2000, new double[]{1, 1, 3}, Intervals.createMinMax( 0,0,0, 512, 512, 200 ), Intervals.createMinMax( 0,0,0, 512, 512, 200 ) );
		sb.addTile( 0, new double[]{0,0,0} );
		sb.addTile( 1, new double[]{10,0,0} );
		ImageJFunctions.show( sb.getImg( 0, 0, 0, 0, 0 ) );
		ImageJFunctions.show( sb.getImg( 0, 0, 0, 1, 0 ) );
	}
}
