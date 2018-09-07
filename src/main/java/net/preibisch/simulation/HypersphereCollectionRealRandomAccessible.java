package net.preibisch.simulation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.gson.GsonBuilder;

import ij.ImageJ;
import net.imglib2.FinalInterval;
import net.imglib2.KDTree;
import net.imglib2.RealInterval;
import net.imglib2.RealLocalizable;
import net.imglib2.RealPoint;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.Sampler;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.neighborsearch.RadiusNeighborSearchOnKDTree;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineRealRandomAccessible;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale;
import net.imglib2.realtransform.Scale3D;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

public class HypersphereCollectionRealRandomAccessible< T extends Type< T > > implements RealRandomAccessible< T >
{

	final private T type;
	final private List< RealLocalizable > sphereCenters;
	final private List< T > sphereValues;
	final private List< Double > sphereRadii;
	final private int numDimensions;
	private KDTree< Integer > kdTree;
	private boolean kdNeedsUpdate;
	private double maxRadius;

	public HypersphereCollectionRealRandomAccessible(int numDimensions, T type)
	{
		this.type = type.copy();
		this.numDimensions = numDimensions;
		sphereCenters = new ArrayList<>();
		sphereRadii = new ArrayList<>();
		sphereValues = new ArrayList<>();
		kdNeedsUpdate = true;
		maxRadius = 0.0;
	}

	public HypersphereCollectionRealRandomAccessible(Collection< RealLocalizable > locations, Collection< Double > radii, Collection< T > values, T type)
	{
		this(locations.iterator().next().numDimensions(), type);
		Iterator< RealLocalizable > itLocations = locations.iterator();
		Iterator< Double > itRadii = radii.iterator();
		Iterator< T > itValues = values.iterator();
		while(itLocations.hasNext())
			addSphere( itLocations.next(), itRadii.next(), itValues.next());
	}

	public void addSphere(RealLocalizable location, Double radius, T value)
	{
		sphereCenters.add( new RealPoint( location ) );
		sphereRadii.add( radius );
		sphereValues.add( value.copy() );
		kdNeedsUpdate = true;
		maxRadius = Math.max( maxRadius, radius );
	}

	private synchronized void prepareKDTree()
	{
		if (!kdNeedsUpdate)
			return;

		// the values of the KDTree nodes are just increasing indices
		// we will later use them to access radii, values, ..
		final AtomicInteger i = new AtomicInteger( 0 );
		kdTree = new KDTree<>( sphereCenters.stream().map( (x) -> i.getAndIncrement() ).collect( Collectors.toList() ), sphereCenters );
		kdNeedsUpdate = false;
	}

	@Override
	public int numDimensions()
	{
		return numDimensions;
	}

	@Override
	public RealRandomAccess< T > realRandomAccess()
	{
		prepareKDTree();
		return new HypersphereCollectionRealRandomAccess();
	}

	@Override
	public RealRandomAccess< T > realRandomAccess(RealInterval interval)
	{
		return realRandomAccess();
	}

	private class HypersphereCollectionRealRandomAccess extends RealPoint implements RealRandomAccess< T >
	{

		final private T type;
		final private RadiusNeighborSearchOnKDTree< Integer > rs;
		
		public HypersphereCollectionRealRandomAccess()
		{
			super(numDimensions);
			this.type = HypersphereCollectionRealRandomAccessible.this.type.copy();
			this.rs = new RadiusNeighborSearchOnKDTree<>( kdTree );
		}

		@Override
		public T get()
		{

			rs.search( this, maxRadius, false );

			if (rs.numNeighbors() < 1)
			{
				type.set( HypersphereCollectionRealRandomAccessible.this.type );
			}
			else
			{
				boolean found = false;
				int minI = Integer.MAX_VALUE;
				for (int i=0; i<rs.numNeighbors(); i++ )
				{
					final Double radius = sphereRadii.get( rs.getSampler( i ).get() );
					final RealLocalizable pos = rs.getPosition( i );
					if (Util.distance( pos, this ) <= radius)
					{
						minI = Math.min( minI, rs.getSampler( i ).get() );
						found = true;
					}
					if (found)
						type.set( sphereValues.get(minI ) );
				}
				if(!found)
					type.set( HypersphereCollectionRealRandomAccessible.this.type );
			}
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
			RealRandomAccess< T > res = new HypersphereCollectionRealRandomAccess();
			res.setPosition( this );
			return res;
		}

	}

}
