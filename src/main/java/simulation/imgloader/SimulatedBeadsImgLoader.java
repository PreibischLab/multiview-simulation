package simulation.imgloader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import mpicbg.spim.data.SpimData;
import mpicbg.spim.data.SpimDataException;
import mpicbg.spim.data.XmlIoSpimData;
import mpicbg.spim.data.registration.ViewRegistration;
import mpicbg.spim.data.registration.ViewRegistrations;
import mpicbg.spim.data.registration.ViewTransform;
import mpicbg.spim.data.registration.ViewTransformAffine;
import mpicbg.spim.data.sequence.Angle;
import mpicbg.spim.data.sequence.Channel;
import mpicbg.spim.data.sequence.FinalVoxelDimensions;
import mpicbg.spim.data.sequence.Illumination;
import mpicbg.spim.data.sequence.ImgLoader;
import mpicbg.spim.data.sequence.MissingViews;
import mpicbg.spim.data.sequence.SequenceDescription;
import mpicbg.spim.data.sequence.TimePoint;
import mpicbg.spim.data.sequence.TimePoints;
import mpicbg.spim.data.sequence.ViewDescription;
import mpicbg.spim.data.sequence.ViewId;
import mpicbg.spim.data.sequence.ViewSetup;
import mpicbg.spim.data.sequence.VoxelDimensions;
import net.imglib2.Cursor;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.type.numeric.real.FloatType;
import simulation.SimulateBeads;

public class SimulatedBeadsImgLoader implements ImgLoader< UnsignedShortType >
{
	final SimulateBeads sb;

	public SimulatedBeadsImgLoader(
			final int[] angles,
			final int axis,
			final int numPoints,
			final Interval rangeSimulation,
			final Interval intervalRender,
			final double[] sigma )
	{
		this( new SimulateBeads( angles, axis, numPoints, rangeSimulation, intervalRender, sigma ) );
	}

	public SimulatedBeadsImgLoader( final SimulateBeads sb ) { this.sb = sb; }

	public SimulateBeads getSimulateBeads() { return sb; }

	@Override
	public RandomAccessibleInterval< UnsignedShortType > getImage( final ViewId view )
	{
		final long[] dim = new long[ sb.getImgs().get( view.getViewSetupId() ).numDimensions() ];

		for ( int d = 0; d < dim.length; ++d )
			dim[ d ] = sb.getImgs().get( view.getViewSetupId() ).dimension( d );

		final Img< UnsignedShortType > img = ArrayImgs.unsignedShorts( dim );

		final Cursor< FloatType > in = sb.getImgs().get( view.getViewSetupId() ).cursor();
		final Cursor< UnsignedShortType > out = img.cursor();

		while ( in.hasNext() )
			out.next().set( Math.round( in.next().get() ) );

		return img;
	}

	@Override
	public UnsignedShortType getImageType() { return new UnsignedShortType(); }

	@Override
	public RandomAccessibleInterval<FloatType> getFloatImage( final ViewId view, boolean normalize )
	{
		if ( normalize )
		{
			return normalize( sb.getImgs().get( view.getViewSetupId() ).copy() );
		}
		else
		{
			return sb.getImgs().get( view.getViewSetupId() ).copy();
		}
	}

	@Override
	public Dimensions getImageSize( final ViewId view ) { return sb.getImgs().get( view.getViewSetupId() ); }

	@Override
	public VoxelDimensions getVoxelSize( final ViewId view ) { return new FinalVoxelDimensions( "pixel", 1, 1, 1 ); }

	protected static final Img< FloatType > normalize( final Img< FloatType > img )
	{
		float min = Float.MAX_VALUE;
		float max = -Float.MAX_VALUE;

		for ( final FloatType t : img )
		{
			final float v = t.get();

			if ( v < min )
				min = v;

			if ( v > max )
				max = v;
		}

		for ( final FloatType t : img )
			t.set( ( t.get() - min ) / ( max - min ) );

		return img;
	}

	public static void save( final SpimData spimData, final String xmlFilename ) throws SpimDataException
	{
		XmlIoSpimData xml = new XmlIoSpimData();
		xml.save( spimData, xmlFilename );
	}

	public static SpimData spimdataExample()
	{
		final int[] angles = new int[]{ 0, 45, 90, 135 };
		final int axis = 0;
		final int numPoints = 1000;
		final double[] sigma = new double[]{ 1, 1, 3 };
		final FinalInterval range = new FinalInterval( 512, 512, 200 );

		return spimdataExample( angles, axis, numPoints, sigma, range );
	}

	public static SpimData spimdataExample( final int[] angles, final int axis, final int numPoints, final double[] sigma, final Interval range )
	{
		final TimePoints timepoints = createTimepoints();
		final ArrayList< ViewSetup > setups = createViewSetups( angles, axis, range );
		final MissingViews missingViews = null;

		final SimulateBeads sb = new SimulateBeads( angles, axis, numPoints, range, range, sigma );

		final SequenceDescription sequenceDescription = new SequenceDescription( timepoints, setups, null, missingViews );
		final ImgLoader< UnsignedShortType > imgLoader = new SimulatedBeadsImgLoader( sb );
		sequenceDescription.setImgLoader( imgLoader );

		// get the minimal resolution of all calibrations
		final double minResolution = 1.0;

		final ViewRegistrations viewRegistrations = createViewRegistrations( sequenceDescription.getViewDescriptions(), minResolution );

		// finally create the SpimData itself based on the sequence description and the view registration
		final SpimData spimData = new SpimData( new File( "" ), sequenceDescription, viewRegistrations );

		return spimData;
	}

	public static ViewRegistrations createViewRegistrations( final Map< ViewId, ViewDescription > viewDescriptionList, final double minResolution )
	{
		final HashMap< ViewId, ViewRegistration > viewRegistrationList = new HashMap< ViewId, ViewRegistration >();
		
		for ( final ViewDescription viewDescription : viewDescriptionList.values() )
			if ( viewDescription.isPresent() )
			{
				final ViewRegistration viewRegistration = new ViewRegistration( viewDescription.getTimePointId(), viewDescription.getViewSetupId() );
				
				final VoxelDimensions voxelSize = viewDescription.getViewSetup().getVoxelSize(); 

				final double calX = voxelSize.dimension( 0 ) / minResolution;
				final double calY = voxelSize.dimension( 1 ) / minResolution;
				final double calZ = voxelSize.dimension( 2 ) / minResolution;
				
				final AffineTransform3D m = new AffineTransform3D();
				m.set( calX, 0.0f, 0.0f, 0.0f, 
					   0.0f, calY, 0.0f, 0.0f,
					   0.0f, 0.0f, calZ, 0.0f );
				final ViewTransform vt = new ViewTransformAffine( "calibration", m );
				viewRegistration.preconcatenateTransform( vt );
				
				viewRegistrationList.put( viewRegistration, viewRegistration );
			}
		
		return new ViewRegistrations( viewRegistrationList );
	}

	public static TimePoints createTimepoints()
	{
		final ArrayList< TimePoint > timepointList = new ArrayList< TimePoint >();
		timepointList.add( new TimePoint( 0 ) );
		return new TimePoints( timepointList );
	}

	public static ArrayList< ViewSetup > createViewSetups( final int[] rotations, final int rotationAxis, final Interval range )
	{
		final ArrayList< Channel > channels = new ArrayList< Channel >();
		channels.add( new Channel( 0, "Channel1" ) );

		final ArrayList< Illumination > illuminations = new ArrayList< Illumination >();
		illuminations.add( new Illumination( 0, "Illum1" ) );

		final ArrayList< Angle > angles = new ArrayList< Angle >();
		for ( int a = 0; a < rotations.length; ++a )
		{
			final Angle angle = new Angle( a, String.valueOf( rotations[ a ] ) );

			final double degrees = rotations[ a ];
			double[] axis = new double[ 3 ];

			axis[ rotationAxis ] = 1;
			angle.setRotation( axis, degrees );

			angles.add( angle );
		}

		final ArrayList< ViewSetup > viewSetups = new ArrayList< ViewSetup >();
		for ( final Channel c : channels )
			for ( final Illumination i : illuminations )
				for ( final Angle a : angles )
				{
					final VoxelDimensions voxelSize = new FinalVoxelDimensions( "pixels", 1, 1, 1 );
					final Dimensions dim = new FinalDimensions( new long[]{ range.max( 0 ) - range.min( 0 ), range.max( 1 ) - range.min( 1 ), range.max( 2 ) - range.min( 2 ) } );
					viewSetups.add( new ViewSetup( viewSetups.size(), null, dim, voxelSize, c, a, i ) );
				}

		return viewSetups;
	}

	public static void main( String[] args ) throws SpimDataException
	{
		final SpimData d = SimulatedBeadsImgLoader.spimdataExample();
		save( d, new File( "simulated.xml" ).getAbsolutePath() );
		System.out.println( "done" );
	}
}
