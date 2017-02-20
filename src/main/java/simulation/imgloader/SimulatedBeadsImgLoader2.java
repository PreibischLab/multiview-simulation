package simulation.imgloader;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import ij.gui.GenericDialog;
import mpicbg.spim.data.SpimData;
import mpicbg.spim.data.legacy.LegacyImgLoaderWrapper;
import mpicbg.spim.data.registration.ViewRegistration;
import mpicbg.spim.data.registration.ViewRegistrations;
import mpicbg.spim.data.registration.ViewTransform;
import mpicbg.spim.data.registration.ViewTransformAffine;
import mpicbg.spim.data.sequence.Angle;
import mpicbg.spim.data.sequence.Channel;
import mpicbg.spim.data.sequence.FinalVoxelDimensions;
import mpicbg.spim.data.sequence.Illumination;
import mpicbg.spim.data.sequence.SequenceDescription;
import mpicbg.spim.data.sequence.Tile;
import mpicbg.spim.data.sequence.TimePoint;
import mpicbg.spim.data.sequence.TimePoints;
import mpicbg.spim.data.sequence.ViewDescription;
import mpicbg.spim.data.sequence.ViewId;
import mpicbg.spim.data.sequence.ViewSetup;
import mpicbg.spim.data.sequence.VoxelDimensions;
import net.imglib2.Dimensions;
import net.imglib2.FinalDimensions;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.Translation;
import net.imglib2.realtransform.Translation3D;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import simulation.RegularTranformHelpers;
import simulation.RegularTranformHelpers.RegularTranslationParameters;
import simulation.SimulateBeads2;

public class SimulatedBeadsImgLoader2 extends LegacyImgLoaderWrapper< UnsignedShortType, LegacySimulatedBeadsImgLoader2 >
{
	
	public SimulatedBeadsImgLoader2(final SimulateBeads2 sb,
			int rotAxis,
			double[] ratations,
			List<double[]> channelShifts,
			List<double[]> illumShifts,
			List<double[]> timeShifts,
			List<double[]> tileShifts
			)
	{
		super( new LegacySimulatedBeadsImgLoader2( sb ) );
		
		for (int tp = 0; tp<timeShifts.size(); tp++)
		{
			sb.addTimepoint( tp, timeShifts.get( tp ) );
			int vsid = 0;
			for (int a = 0; a<ratations.length; a++)
				for (int channel = 0; channel<channelShifts.size(); channel++)
					for (int illum = 0; illum<illumShifts.size(); illum++)
						for (int tile = 0; tile<tileShifts.size(); tile++)
						{
							
							sb.addAngle( a, rotAxis, ratations[a] );
							sb.addChannel( channel, channelShifts.get( channel ) );
							sb.addIllumination( illum, illumShifts.get( illum ) );
							sb.addTile( tile, tileShifts.get( tile ) );
							addViewId( vsid, tp, a, channel, tile, illum );
							vsid++;
						}
		}
	}

	public SimulatedBeadsImgLoader2( final SimulateBeads2 sb)
	{
		super( new LegacySimulatedBeadsImgLoader2( sb ) );
	}
	
	public SimulateBeads2 getSimulateBeads() { return legacyImgLoader.getSimulateBeads(); }
	
	public void addViewId(int idx, int tp, int angle, int channel, int tile, int illumination)
	{
		legacyImgLoader.addViewId( idx, tp, angle, channel, tile, illumination );
	}
	
	public static ArrayList< ViewSetup > createViewSetupsFromImgLoader(SimulatedBeadsImgLoader2 imgLoader)
	{
		final ArrayList< ViewSetup > viewSetups = new ArrayList< ViewSetup >();
		SimulateBeads2 sb = imgLoader.getSimulateBeads();
		int vid = 0;
			for (int a : new ArrayList<>(sb.angleTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
				for (int channel : new ArrayList<>(sb.channelTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
					for (int illum : new ArrayList<>(sb.illumTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
						for (int tile : new ArrayList<>(sb.tileTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
						{
							final VoxelDimensions voxelSize = new FinalVoxelDimensions( "pixels", 1, 1, 1 );
							final Dimensions dim = imgLoader.getSimulateBeads().getImg( 0, a, channel, tile, illum );
							
							viewSetups.add( new ViewSetup( vid++, Integer.toString( vid ), dim, voxelSize, new Tile(tile), new Channel(channel), new Angle(a), new Illumination(illum) ) );
						}
		
		return viewSetups;
	}
		
	public static TimePoints createTimePointsFromImgLoader(SimulatedBeadsImgLoader2 imgLoader)
	{
		final ArrayList< TimePoint > timepointList = new ArrayList< TimePoint >();
		SimulateBeads2 sb = imgLoader.getSimulateBeads();
		for (int tp : new ArrayList<>(sb.tpTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
			timepointList.add( new TimePoint( tp ) );
		return new TimePoints( timepointList );
			
	}
	
	public static ViewRegistrations createViewRegistrationsFromImgLoader(SimulatedBeadsImgLoader2 imgLoader, double relativeTileError, boolean centerAngles)
	{
		SimulateBeads2 sb = imgLoader.getSimulateBeads();
		final HashMap< ViewId, ViewRegistration > viewRegistrationList = new HashMap< ViewId, ViewRegistration >();
		final RealInterval tilesExtent = sb.getTilesExtent();
		double[] center = new double[tilesExtent.numDimensions()];
		tilesExtent.realMax( center );
		
		for (int d = 0; d < center.length; d++)
			center[d] /= -2.0;
		
		for (int tp : new ArrayList<>(sb.tpTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
		{
			int vid = 0;
			for (int a : new ArrayList<>(sb.angleTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
			for (int channel : new ArrayList<>(sb.channelTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
				for (int illum : new ArrayList<>(sb.illumTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
					for (int tile : new ArrayList<>(sb.tileTransforms.keySet()).stream().sorted().collect( Collectors.toList() ))
					{
						final VoxelDimensions voxelSize = new FinalVoxelDimensions( "pixels", 1, 1, 1 );
						final Dimensions dim = imgLoader.getSimulateBeads().getImg( 0, a, channel, tile, illum );
						
						
						
						
						final ViewRegistration viewRegistration = new ViewRegistration( tp, vid++ );
						

						final double calX = voxelSize.dimension( 0 );
						final double calY = voxelSize.dimension( 1 );
						final double calZ = voxelSize.dimension( 2 );
						
						final AffineTransform3D m = new AffineTransform3D();
						m.set( calX, 0.0f, 0.0f, 0.0f, 
							   0.0f, calY, 0.0f, 0.0f,
							   0.0f, 0.0f, calZ, 0.0f );
						final ViewTransform vt = new ViewTransformAffine( "calibration", m );
						viewRegistration.preconcatenateTransform( vt );
						
						double[] translation = sb.tileTransforms.get( tile ).inverse().getTranslation();
						for (int i = 0; i < translation.length; i++)
							translation[i] *= relativeTileError;
						
						AffineTransform3D atr = new AffineTransform3D();
						atr.translate( translation );
						final ViewTransform vtT = new ViewTransformAffine("translation" , atr.copy() );
						viewRegistration.preconcatenateTransform( vtT );
						
						if (centerAngles)
						{
							AffineTransform3D centerTr = new AffineTransform3D();
							centerTr.translate( center );
							final ViewTransform centerVt = new ViewTransformAffine("center angle" , centerTr.copy() );
							viewRegistration.preconcatenateTransform( centerVt );
						}
						
						final ViewTransform vtA = new ViewTransformAffine("rotation" , sb.angleTransforms.get( a ));
						viewRegistration.preconcatenateTransform( vtA );
						
						
						
						
						
						viewRegistrationList.put( viewRegistration, viewRegistration );
						
							
					}
		}
		return new ViewRegistrations( viewRegistrationList );
	}
	
	public static SpimData createSpimData(
			int numPoints,
			double[] sigma,
			Interval rangeSimulation,
			Interval intervalRender,
			int rotAxis,
			double[] ratations,
			List<double[]> channelShifts,
			List<double[]> illumShifts,
			List<double[]> timeShifts,
			List<double[]> tileShifts,
			double relativeTileError,
			boolean centerAngles
	)
	{
		SimulateBeads2 sb = new SimulateBeads2( numPoints, sigma, rangeSimulation, intervalRender );
		SimulatedBeadsImgLoader2 loader = new SimulatedBeadsImgLoader2( sb, rotAxis, ratations, channelShifts, illumShifts, timeShifts, tileShifts );
		
		
		SequenceDescription sd = new SequenceDescription( createTimePointsFromImgLoader( loader ), createViewSetupsFromImgLoader( loader ) );
		
		sd.setImgLoader( loader );
		SpimData res = new SpimData( new File(""), sd, createViewRegistrationsFromImgLoader( loader, relativeTileError, centerAngles ) );
		
		return res;
		
	}
	
	
	public static SpimData createSpimDataFromUserInput()
	{
		GenericDialog gd = new GenericDialog( "Simulated SpimData Parameters" );
		
		gd.addNumericField( "number of points", 2000, 0 );
		gd.addStringField( "Gaussian sigmas", "1,1,3" );
		
		gd.addMessage( "Please enter Intervals in the form: min1,...,mini,max1,...maxi" );
		gd.addStringField( "Simulation Interval", "-512,-512,-512,512,512,512", 35 );
		gd.addStringField( "Viewport Interval", "0,0,0,256,256,100", 35 );
		
		
		gd.addNumericField( "Number of TimePoints", 1, 0 );
		gd.addNumericField( "Number of Channels", 1, 0 );
		gd.addNumericField( "Number of Illuminations", 1, 0 );
		
		gd.addNumericField( "Number of Tiles X", 1, 0 );
		gd.addNumericField( "Number of Tiles Y", 1, 0 );
		gd.addNumericField( "Tile Overlap", 0.15, 3);
		gd.addMessage( "Actual tile shifts will be multiplied by relative error in metadata" );
		gd.addNumericField( "Tile Metadata Error", 1, 3);
		
		gd.addStringField( "Angles", "0,90" );
		gd.addNumericField( "Angle Axis of Rotation", 1, 0 );
		
		gd.addCheckbox( "center angles before rotation", true );
		
		
		gd.showDialog();
		if (gd.wasCanceled())
			return null;
		
		int numPoints = (int) gd.getNextNumber();
		String[] sigmasString = gd.getNextString().split( "," );
		double[] sigmas = new double[sigmasString.length];
		{
			AtomicInteger i = new AtomicInteger( 0 );
			Arrays.asList( sigmasString ).forEach( (s) -> sigmas[i.getAndIncrement()] = Double.parseDouble( s ) );
		}
		String[] simMinMaxString = gd.getNextString().split( "," );
		long[] simMinMax = new long[simMinMaxString.length];
		{
			AtomicInteger i = new AtomicInteger( 0 );
			Arrays.asList( simMinMaxString ).forEach( (s) -> simMinMax[i.getAndIncrement()] = Long.parseLong( s ) );
		}
		String[] vpMinMaxString = gd.getNextString().split( "," );
		long[] vpMinMax = new long[vpMinMaxString.length];
		{
			AtomicInteger i = new AtomicInteger( 0 );
			Arrays.asList( vpMinMaxString ).forEach( (s) -> vpMinMax[i.getAndIncrement()] = Long.parseLong( s ) );
		}
		
		int numTP = Math.max(1, (int) gd.getNextNumber());
		int numChannels = Math.max(1, (int) gd.getNextNumber());
		int numIllums = Math.max(1, (int) gd.getNextNumber());
		
		int numTilesX = Math.max(1, (int) gd.getNextNumber());
		int numTilesY = Math.max(1, (int) gd.getNextNumber());
		
		double overlap =  gd.getNextNumber();
		double relError =  gd.getNextNumber();
		
		String[] anglesString = gd.getNextString().split( "," );
		double[] angles = new double[anglesString.length];
		{
			AtomicInteger i = new AtomicInteger( 0 );
			Arrays.asList( anglesString ).forEach( (s) -> angles[i.getAndIncrement()] = Double.parseDouble( s ) );
		}
		
		int rotAxis = (int) gd.getNextNumber();
		
		boolean center = gd.getNextBoolean();
		
		
		
		GenericDialog gd2 = new GenericDialog("View Errors");
		
		for (int i = 0 ; i < numTP; i++)
			gd2.addStringField( "TP " + i + " error", "0,0,0" );
		
		for (int i = 0 ; i < numChannels; i++)
			gd2.addStringField( "Channel " + i + " error", "0,0,0" );
		
		for (int i = 0 ; i < numIllums; i++)
			gd2.addStringField( "Illumination  " + i + " error", "0,0,0" );
		
		gd2.showDialog();
		if (gd2.wasCanceled())
			return null;
		
		List<double[]> channelShifts = new ArrayList<>();
		List<double[]> illumShifts = new ArrayList<>();
		List<double[]> timeShifts = new ArrayList<>();
		
		for (int i = 0 ; i < numTP; i++)
		{
			String[] shiftString = gd2.getNextString().split( "," );
			double[] shift = new double[shiftString.length];
			{
				AtomicInteger j = new AtomicInteger( 0 );
				Arrays.asList( shiftString ).forEach( (s) -> shift[j.getAndIncrement()] = Double.parseDouble( s ) );
			}
			timeShifts.add( shift );
		}
		
		for (int i = 0 ; i < numChannels; i++)
		{
			String[] shiftString = gd2.getNextString().split( "," );
			double[] shift = new double[shiftString.length];
			{
				AtomicInteger j = new AtomicInteger( 0 );
				Arrays.asList( shiftString ).forEach( (s) -> shift[j.getAndIncrement()] = Double.parseDouble( s ) );
			}
			channelShifts.add( shift );
		}
		
		for (int i = 0 ; i < numIllums; i++)
		{
			String[] shiftString = gd2.getNextString().split( "," );
			double[] shift = new double[shiftString.length];
			{
				AtomicInteger j = new AtomicInteger( 0 );
				Arrays.asList( shiftString ).forEach( (s) -> shift[j.getAndIncrement()] = Double.parseDouble( s ) );
			}
			illumShifts.add( shift );
		}
		
		RegularTranslationParameters params = new RegularTranslationParameters();
		params.nDimensions = 3;
		params.alternating = new boolean[] {true, true, true};
		params.dimensionOrder = new int[] {0, 1, 2};
		params.increasing = new boolean[] {true, true, true};
		params.overlaps = new double[] {overlap, overlap, overlap};
		params.nSteps = new int[] {numTilesX, numTilesY, 1};
		List< AffineTransform3D > generateRegularGrid = RegularTranformHelpers.generateRegularGrid( params, Intervals.createMinMax( vpMinMax ) );
		
		List<double[]> tileShifts = new ArrayList<>();
		generateRegularGrid.forEach( ( t ) -> tileShifts.add( t.getTranslation() ));
	
		
		return createSpimData( numPoints, sigmas, Intervals.createMinMax( simMinMax ), Intervals.createMinMax( vpMinMax ), rotAxis, angles, channelShifts, illumShifts, timeShifts, tileShifts, relError, center );
	}
	
	public static void main(String[] args)
	{
		SpimData sd2 = createSpimDataFromUserInput();
		RandomAccessibleInterval< UnsignedShortType > img2 = (RandomAccessibleInterval< UnsignedShortType >) sd2.getSequenceDescription().getImgLoader().getSetupImgLoader( 0 ).getImage( 0, null );
		ImageJFunctions.show( img2 );
		if (true)
			return;
		
		List<double[]> channelShifts = new ArrayList<>();
		channelShifts.add( new double[3] );
		List<double[]> illumShifts = new ArrayList<>();
		illumShifts.add( new double[3] );
		List<double[]> timeShifts = new ArrayList<>();
		timeShifts.add( new double[3] );
		timeShifts.add( new double[3] );
		List<double[]> tileShifts = new ArrayList<>();
		tileShifts.add( new double[3] );
		
		SpimData sd = createSpimData( 1000, new double[] {1, 1, 3 }, Intervals.createMinMax( 0,0,0,256,256,100 ),
				Intervals.createMinMax( 0,0,0,256,256,100 ), 1, new double[] {0, 90}, channelShifts, illumShifts, timeShifts, tileShifts , 0.9, true);
		
		RandomAccessibleInterval< UnsignedShortType > img = (RandomAccessibleInterval< UnsignedShortType >) sd.getSequenceDescription().getImgLoader().getSetupImgLoader( 0 ).getImage( 1, null );
		ImageJFunctions.show( img );
	}

}
