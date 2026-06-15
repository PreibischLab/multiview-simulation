/*-
 * #%L
 * Library for simulating a multi-view acquisition including
 * attenuation, convolution, reduced sampling and poission noise.
 * %%
 * Copyright (C) 2014 - 2026 Multiview Simulation developers.
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
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Executors;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class PerlinNoiseRealRandomAccessible<T extends RealType< T >> implements RealRandomAccessible< T >
{

	private final int numDim;
	private final double[] scales;
	private final int[] loopExtents;
	private final int nVectors;
	private final T type;

	private final List<double[]> gradients;
	private final List<Integer> permuatation;

	public PerlinNoiseRealRandomAccessible(T type, double[] scales, int[] loopExtents, int nVectors, Random rng)
	{
		this.type = type.copy();
		this.nVectors = nVectors;

		this.scales = scales.clone();
		this.loopExtents = loopExtents.clone();
		this.numDim = scales.length;

		gradients = new ArrayList<>();
		permuatation = new ArrayList<>();
		for (int i=0; i<nVectors; i++)
		{
			gradients.add( randomUnitVector( numDim, rng ) );
			permuatation.add( i );
		}
		Collections.shuffle( permuatation, rng );
	}

	@Override
	public int numDimensions()
	{
		return numDim;
	}

	@Override
	public T getType()
	{
		return type;
	}

	@Override
	public RealRandomAccess< T > realRandomAccess()
	{
		return new PerlinNoiseRealRandomAccess();
	}

	@Override
	public RealRandomAccess< T > realRandomAccess(RealInterval interval)
	{
		return realRandomAccess();
	}
	
	public static double[] randomUnitVector(int nDim, Random rng)
	{
		double[] res = new double[nDim];
		double sSum = 0.0;
		for (int d=0; d<nDim; d++)
		{
			res[d] = rng.nextGaussian();
			sSum += res[d] * res[d];
		}
		for (int d=0; d<nDim; d++)
			res[d] /= Math.sqrt( sSum );
		return res;
	}
	
	private class PerlinNoiseRealRandomAccess extends RealPoint implements RealRandomAccess< T >
	{

		private final T type;
		private final List<int[]> neighbors;

		public PerlinNoiseRealRandomAccess()
		{
			super(numDim);
			this.type = PerlinNoiseRealRandomAccessible.this.type.copy();
			neighbors = neighborOffsets( numDim );
		}
		
		@Override
		public T get()
		{
			double[] posInGrid = position.clone();
			int[] posInGridInt = new int[posInGrid.length];
			for (int d=0; d<numDim; d++)
			{
				posInGrid[d] = fFloorMod( (float) ( posInGrid[d] / scales[d] ), loopExtents[d]);
				posInGridInt[d] = (int) Math.floor( posInGrid[d] );
			}

			double[] offs = null;
			List<Double> dots = new ArrayList<>();
			for(int[] off : neighbors)
			{
				int[] neighborPos = posInGridInt.clone();
				double[] dist = posInGrid.clone();
				for (int d=0; d<numDim; d++)
				{
					neighborPos[d] = neighborPos[d] + off[d];
					dist[d] = posInGrid[d] - neighborPos[d];
					neighborPos[d] = neighborPos[d] % loopExtents[d];
				}
				
				if (offs == null)
					offs = dist;
				int idx = flatIndex( neighborPos, loopExtents ) % nVectors;
				double[] grad = gradients.get( permuatation.get( idx ) );
				
				dots.add(dot(grad, dist));
			}

			type.setReal( interpolateSmoothstep( dots, offs ) );
			return type;
		}

		@Override
		public T getType()
		{
			return type;
		}

		@Override
		public RealRandomAccess< T > copy()
		{
			return copyRealRandomAccess();
		}

		@Override
		public RealRandomAccess< T > copyRealRandomAccess()
		{
			PerlinNoiseRealRandomAccess ra = new PerlinNoiseRealRandomAccess();
			ra.setPosition( this );
			return ra;
		}
		
	}

	private static double interpolateSmoothstep(List<Double> dots, double[] offs)
	{
		List<Double> interpolated = new ArrayList<>();
		interpolated.addAll( dots );
		for (int d=offs.length-1; d>=0; d--)
		{
			List<Double> interpolatedTemp = new ArrayList<>();
			Iterator< Double > it = interpolated.iterator();
			for (int i=0; i<(int)Math.pow( 2, d );i++)
			{
				Double a1 = it.next();
				Double a2 = it.next();
				interpolatedTemp.add( smoothstep(a1, a2, offs[d]) );
			}
			interpolated = interpolatedTemp;
		}
		return interpolated.get( 0 );
	}
	
	private static double smoothstep(double a1, double a2, double p)
	{
		double sstep = Math.pow( p, 3 ) * (10 - 15 * p + 6* Math.pow( p, 2 ));
		sstep = Math.min( 1, Math.max( sstep, 0 ) );
		return (1.0 - sstep) * a1 + sstep * a2;
	}
	
	private static double dot(double[] a, double[] b)
	{
		double dot = 0;
		for (int d=0; d<a.length; d++)
			dot += a[d] * b[d];
		return dot;
	}

	private static int flatIndex(int[] pos, int[] dim)
	{
		int cp = 1;
		int idx = 0;
		for (int d=0; d<pos.length; d++)
		{
			idx += pos[d] * cp;
			cp += cp * dim[d];
		}
		return idx;
	}

	private static List<int[]> neighborOffsets(int numDim)
	{
		List<int[]> res = new ArrayList<>();
		for (int i=0; i<Math.pow(2, numDim); i++)
		{
			int it = i;
			int[] off = new int[numDim];
			for (int d=numDim-1;d>=0;d--)
			{
				off[numDim - 1 - d] = it / (int) Math.pow( 2, d );
				it = it % (int) Math.pow( 2, d );
			}
			res.add( off );
		}
		return res;
	}
	
	
	public static double fFloorMod(float a, float b)
	{
		float mod = a % b;
		return mod < 0 ? mod + b : mod;
	}

	public static void main(String[] args)
	{
		Random rng = new Random(42);
		PerlinNoiseRealRandomAccessible< FloatType > rrable = new PerlinNoiseRealRandomAccessible<>( new FloatType(), new double[] {512,1500,256}, new int[] {15, 15, 15}, 100, rng );
	
		IntervalView< FloatType > view = Views.interval( Views.raster( rrable ), new FinalInterval( new long[] {2048, 2048, 256} ) );
		new ImageJ();
		ImageJFunctions.show( view , Executors.newFixedThreadPool( Runtime.getRuntime().availableProcessors() ));
	}

}
