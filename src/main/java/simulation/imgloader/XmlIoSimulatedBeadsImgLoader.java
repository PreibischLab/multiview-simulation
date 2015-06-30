package simulation.imgloader;

import static mpicbg.spim.data.XmlKeys.IMGLOADER_FORMAT_ATTRIBUTE_NAME;

import java.io.File;

import mpicbg.spim.data.XmlHelpers;
import mpicbg.spim.data.generic.sequence.AbstractSequenceDescription;
import mpicbg.spim.data.generic.sequence.ImgLoaderIo;
import mpicbg.spim.data.generic.sequence.XmlIoBasicImgLoader;
import net.imglib2.FinalInterval;

import org.jdom2.Element;

import simulation.SimulateBeads;

@ImgLoaderIo( format = "spimreconstruction.simulatedbeads", type = SimulatedBeadsImgLoader.class )
public class XmlIoSimulatedBeadsImgLoader implements XmlIoBasicImgLoader< SimulatedBeadsImgLoader >
{
	public static final String DIRECTORY_TAG = "imagedirectory";
	public static final String MASTER_FILE_TAG = "masterfile";
	public static final String IMGLIB2CONTAINER_PATTERN_TAG = "imglib2container";

	@Override
	public Element toXml( final SimulatedBeadsImgLoader imgLoader, final File basePath )
	{
		final Element elem = new Element( "ImageLoader" );
		elem.setAttribute( IMGLOADER_FORMAT_ATTRIBUTE_NAME, this.getClass().getAnnotation( ImgLoaderIo.class ).format() );

		final SimulateBeads sb = imgLoader.getSimulateBeads();

		final int[] rangeSimulation = new int[ 6 ];
		final int[] intervalRender = new int[ 6 ];

		for ( int d = 0; d < 3; ++d )
		{
			rangeSimulation[ d ] = (int)sb.rangeSimulation.min( d );
			intervalRender[ d ] = (int)sb.intervalRender.min( d );

			rangeSimulation[ d + 3 ] = (int)sb.rangeSimulation.max( d );
			intervalRender[ d + 3 ] = (int)sb.intervalRender.max( d );
		}

		elem.addContent( XmlHelpers.intArrayElement( "rotation_angles", sb.angles ) );
		elem.addContent( XmlHelpers.intElement( "axis", sb.axis ) );
		elem.addContent( XmlHelpers.intElement( "num_points", sb.numPoints ) );
		elem.addContent( XmlHelpers.intArrayElement( "range_simulation", rangeSimulation ) );
		elem.addContent( XmlHelpers.intArrayElement( "interval_render", intervalRender ) );
		elem.addContent( XmlHelpers.doubleArrayElement( "sigma", sb.sigma ) );

		return elem;
	}

	@Override
	public SimulatedBeadsImgLoader fromXml(
			final Element elem, File basePath,
			final AbstractSequenceDescription<?, ?, ?> sequenceDescription )
	{
		final int[] angles = XmlHelpers.getIntArray( elem, "rotation_angles" );
		final int axis = XmlHelpers.getInt( elem, "axis" );
		final int numPoints = XmlHelpers.getInt( elem, "num_points" );
		final int[] rs = XmlHelpers.getIntArray( elem, "range_simulation" );
		final int[] ir = XmlHelpers.getIntArray( elem, "interval_render" );
		final double[] sigma = XmlHelpers.getDoubleArray( elem, "sigma" );

		final FinalInterval rangeSimulation = new FinalInterval( new long[]{ rs[ 0 ], rs[ 1 ], rs[ 2 ] }, new long[]{ rs[ 3 ], rs[ 4 ], rs[ 5 ] } );
		final FinalInterval intervalRender = new FinalInterval( new long[]{ ir[ 0 ], ir[ 1 ], ir[ 2 ] }, new long[]{ ir[ 3 ], ir[ 4 ], ir[ 5 ] } );

		final SimulateBeads sb = new SimulateBeads( angles, axis, numPoints, rangeSimulation, intervalRender, sigma );
		return new SimulatedBeadsImgLoader( sb );
	}
}
