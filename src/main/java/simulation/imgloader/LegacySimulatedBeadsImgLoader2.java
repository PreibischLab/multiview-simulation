package simulation.imgloader;

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
import simulation.SimulateBeads;
import simulation.SimulateBeads2;

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
