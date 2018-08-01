package net.preibisch.simulation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import net.imglib2.RandomAccessible;
import net.imglib2.RealInterval;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.Sampler;
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
		public Sampler< T > copy()
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
