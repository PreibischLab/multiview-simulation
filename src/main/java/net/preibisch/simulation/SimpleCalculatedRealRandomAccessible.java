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
import java.util.Collection;
import java.util.List;

import net.imglib2.RandomAccessible;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.type.Type;

public class SimpleCalculatedRealRandomAccessible<T extends Type< T >> implements RealRandomAccessible< T >
{
	
	final private T type;
	final private List<RealRandomAccessible<T>> in;
	final private Calculation< T > calc;
	final private int numDimensions;
	
	@FunctionalInterface
	interface Calculation<T>
	{
		void calculate(T out, Collection<T> in);
	}

	public SimpleCalculatedRealRandomAccessible(T type, Calculation<T> calc, RealRandomAccessible< T > ... in)
	{
		this.type = type.copy();
		this.in = new ArrayList<>();
		for (RealRandomAccessible< T > i : in)
			this.in.add( i );
		this.numDimensions = in[0].numDimensions();
		this.calc = calc;
	}

	@Override
	public int numDimensions()
	{
		return numDimensions;
	}

	@Override
	public T getType()
	{
		return type;
	}

	@Override
	public RealRandomAccess< T > realRandomAccess()
	{
		return new CalculatedRealRandomAccess();
	}

	@Override
	public RealRandomAccess< T > realRandomAccess(RealInterval interval)
	{
		return realRandomAccess();
	}

	private class CalculatedRealRandomAccess extends RealPoint implements RealRandomAccess< T >
	{

		final private T type;
		final private List<T> inTypes;
		final private List<RealRandomAccess< T >> inRAs;

		public CalculatedRealRandomAccess()
		{
			super(numDimensions);
			this.type = SimpleCalculatedRealRandomAccessible.this.type.copy();
			inTypes = new ArrayList<>();
			inRAs = new ArrayList<>();
			for (int i=0; i<in.size(); i++)
			{
				inRAs.add( in.get( i ).realRandomAccess() );
				inTypes.add( inRAs.get( i ).get().copy() );
			}
		}

		@Override
		public T get()
		{
			for (int i=0; i<inRAs.size(); i++)
			{
				inRAs.get( i ).setPosition( this );
				inTypes.get( i ).set( inRAs.get( i ).get() );
			}
			calc.calculate( type, inTypes );
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
			return new CalculatedRealRandomAccess();
		}
		
	}
	

}
