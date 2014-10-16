//
//  main.m
//  QGDampedSlab
//
//  Created by Jeffrey J. Early on 10/17/12.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLOceanKit/GLOceanKit.h>

int main (int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		GLFloat H1 = 50;
		GLFloat H2 = 800;
		GLFloat dRho1 = 1E-3;
		GLFloat dRho2 = 1E-3;
		GLFloat latitude = 24;
		
        GLFloat domainWidth = 100e3; // m
        NSUInteger nPoints = 256;
        NSUInteger aspectRatio = 1;
        NSUInteger floatFraction = 2;
        
		GLFloat maxTime = 40*86400;
		GLFloat dampingOrder=2; // order of the damping operator. Order 1 is harmonic, order 2 is biharmonic, etc.
		GLFloat dampingTime=3600; // e-folding time scale of the Nyquist frequency.
		GLFloat linearDampingTime = 4*86400;
		
		// Standard constants
		GLFloat f0 = 2 * 7.2921E-5 * sin( latitude*M_PI/180. );
		GLFloat g = 9.81;
		GLFloat R = 6.371e6;
		GLFloat beta = 2 * 7.2921E-5 * cos( latitude*M_PI/180. ) / R;
		GLFloat gprime = dRho1 * g;
		GLFloat rho_water = 1025;
        GLFloat h2 = dRho2*H2;
		GLFloat L_1 = sqrt( gprime * H1) / f0;
		GLFloat L_2 = sqrt( gprime * H2) / f0;

		
        
        GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:nPoints domainMin:-domainWidth/2.0 length:domainWidth];
        xDim.name = @"x"; xDim.units = @"m";
        GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:nPoints/aspectRatio domainMin:-domainWidth/(2.0*aspectRatio) length: domainWidth/aspectRatio];
        yDim.name = @"y"; yDim.units = @"m";
        GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
        tDim.name = @"time"; tDim.units = @"s";
        GLEquation *equation = [[GLEquation alloc] init];
        
        /************************************************************************************************/
        /*		Spin up the lower (QG) layer															*/
        /************************************************************************************************/
        
        NSURL *restartFile = [[NSURL fileURLWithPath: [NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) firstObject]] URLByAppendingPathComponent:@"QGTurbulence.nc"];
        NSFileManager *fileManager = [[NSFileManager alloc] init];
        
        if (![fileManager fileExistsAtPath: restartFile.path])
        {
            Quasigeostrophy2D *qg = [[Quasigeostrophy2D alloc] initWithDimensions: @[xDim, yDim] depth: h2 latitude: latitude equation: equation];
            
            qg.shouldUseBeta = NO;
            qg.shouldUseSVV = YES;
            qg.shouldAntiAlias = NO;
            qg.shouldForce = YES;
            qg.forcingFraction = 16;
            qg.forcingWidth = 1;
            qg.f_zeta = 10;
            qg.forcingDecorrelationTime = HUGE_VAL;
            qg.thermalDampingFraction = 0.0;
            qg.frictionalDampingFraction = 2.0;
            
            qg.outputFile = restartFile;
            qg.shouldAdvectFloats = NO;
            qg.shouldAdvectTracer = NO;
            qg.outputInterval = 10*86400.;
            
            [qg runSimulationToTime: 700*86400];
        }
        
        /************************************************************************************************/
        /*		Create the integrator for the unforced QG layer											*/
        /************************************************************************************************/
        
        Quasigeostrophy2D *qg = [[Quasigeostrophy2D alloc] initWithFile:restartFile resolutionDoubling:NO equation: equation];
        qg.shouldForce = NO;
        
        /************************************************************************************************/
        /*		Create the initial conditions for the slab layer                                        */
        /************************************************************************************************/
        
        NSArray *spatialDimensions = @[xDim, yDim];
        GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
        GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
        
        GLFunction *h0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
        GLFunction *u0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
        GLFunction *v0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
        
        h0 = [h0 setValue: H1 atIndices: @":,:"];
        [u0 zero];
        [v0 zero];
		
        /************************************************************************************************/
        /*		Read the winds from file																*/
        /************************************************************************************************/
        
        // If all goes well, the variable t will be identified as the coordinated variable and therefore turned into a dimension, leaving only u and v.
        GLNetCDFFile *winds = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Documents/Models/QGDampedSlab/winds.nc"] forEquation: equation];
        GLFunction *u_wind = winds.variables[0];
        GLFunction *v_wind = winds.variables[1];
        
#warning For some reason I get total crap if I don't convert from a NetCDFVariable to a regular variable.
        u_wind = [GLFunction functionFromFunction: u_wind];
        v_wind = [GLFunction functionFromFunction: v_wind];
        
        GLFloat rho_air = 1.25;
        GLFloat tau_scale = 1./( rho_water*H1*beta*beta*L_2*L_2*L_2); // Nondimensionalization factor
        GLFunction *speed_wind = [[[u_wind times: u_wind] plus: [v_wind times: v_wind]] sqrt];
        GLFunction *dragCoefficient = [[[speed_wind scalarMultiply: 0.071] scalarAdd: 0.50] scalarMultiply: 1e-3];
        GLFunction *tau_x = [[[[u_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air] times: @(tau_scale)];
        GLFunction *tau_y = [[[[v_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air] times: @(tau_scale)];
        
        [tau_x solve];
        [tau_y solve];
		
		/************************************************************************************************/
		/*		Plop down floats                                                                        */
		/************************************************************************************************/
		
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: xDim.gridType nPoints: xDim.nPoints/floatFraction domainMin: xDim.domainMin length:xDim.domainLength];
        xFloatDim.name = @"x-float";
        GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: yDim.gridType nPoints: yDim.nPoints/floatFraction domainMin: yDim.domainMin length:yDim.domainLength];
        yFloatDim.name = @"y-float";
        
        NSArray *floatDims = @[xFloatDim, yFloatDim];
        GLFunction *xPosition = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDims forEquation: equation];
        GLFunction *yPosition = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDims forEquation: equation];
		
        /************************************************************************************************/
        /*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
        /************************************************************************************************/
        
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Volumes/Data/QGPlusSlab/WindForcedFPlane2.nc"] forEquation: equation overwriteExisting: YES];
        
        GLMutableVariable *hHistory = [h0 variableByAddingDimension: tDim];
        hHistory.name = @"h";
        hHistory.units = @"m";
        hHistory = [netcdfFile addVariable: hHistory];
        
        GLMutableVariable *uHistory = [u0 variableByAddingDimension: tDim];
        uHistory.name = @"u";
        uHistory.units = @"m/s";
        uHistory = [netcdfFile addVariable: uHistory];
        
        GLMutableVariable *vHistory = [v0 variableByAddingDimension: tDim];
        vHistory.name = @"v";
        vHistory.units = @"m/s";
        vHistory = [netcdfFile addVariable: vHistory];
		
//		GLMutableVariable *sshHistory = [ssh variableByAddingDimension: tDim];
//		sshHistory.name = @"SSH";
//		sshHistory = [netcdfFile addVariable: sshHistory];
//		
//		GLMutableVariable *sshFDHistory, *rvHistory;
//		if (writeExtrasToFile) {
//			sshFDHistory = [sshFD variableByAddingDimension: tDim];
//			sshFDHistory.name = @"SSH_FD";
//			sshFDHistory = [netcdfFile addVariable: sshFDHistory];
//			
//			rvHistory = [rv variableByAddingDimension: tDim];
//			rvHistory.name = @"RV";
//			rvHistory = [netcdfFile addVariable: rvHistory];
//		}
		
		GLMutableVariable *xPositionHistory = [xPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
		
		GLMutableVariable *yPositionHistory = [yPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
        /************************************************************************************************/
        /*		Determine an appropriate time step based on the CFL condition.							*/
        /************************************************************************************************/
        
        CGFloat cfl = 0.5;
        GLFloat timeStep = cfl * xDim.sampleInterval / sqrt(g*H1);
        timeStep = 300;
		
		/************************************************************************************************/
		/*		Determine an appropriate time step based on the CFL condition.							*/
		/************************************************************************************************/
        
        GLFloat u_scale = f0/(beta*L_2);
        
        GLFloat div_scale = u_scale*(L_1*L_1)/(L_2*L_2);
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[u0, v0, h0] stepSize: timeStep fFromTY:^(GLScalar *t, NSArray *yNew) {
            GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[t]];
            GLFunction *tx = interp.result[0];
            GLFunction *ty = interp.result[1];
            
            GLFunction *u = yNew[0];
            GLFunction *v = yNew[1];
            GLFunction *h = yNew[2];
            
            GLFunction *fu = [[[v times: @(u_scale)] plus: tx] minus: [[[[h x] scalarMultiply: g1] minus: [u differentiateWithOperator:dampingOp]] spatialDomain]];
            GLFunction *fv = [[[u times: @(-f0)] plus: [[ty times: @(1/rho_water)] dividedBy: h]] minus: [[[[h y] scalarMultiply: g1] minus: [v differentiateWithOperator:dampingOp]] spatialDomain]];
            GLFunction *fh = [[[[u x] plus: [v y]] times: h] negate];
            
            NSArray *f = @[fu, fv, fh];
            return f;
        }];
		
		
		srand(1);
		for (GLFloat time = 1/T_QG; time < maxTime/T_QG; time += 1/T_QG)
		{
			@autoreleasepool {
				yin = [integrator stepForward: yin toTime: time];
				
				zeta = yin[0];
				xPosition = yin[1];
				yPosition = yin[2];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime*T_QG), integrator.lastStepSize*T_QG);
				// We're using spectral code, so it's possible (and is in fact the case) that the variable is not in the spatial domain.
				[tDim addPoint: @(integrator.currentTime)];
				
				GLVariable *etaFD = [zeta diff: @"inverseLaplacianMinusOne"];
				GLVariable *eta = [etaFD spatialDomain];
				GLVariable *rv2 = [[etaFD diff: @"harmonicOperator"] spatialDomain];
				
				[sshHistory concatenateWithLowerDimensionalVariable: eta alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[netcdfFile waitUntilAllOperationsAreFinished];
				
				if (writeExtrasToFile) {
					[sshFDHistory concatenateWithLowerDimensionalVariable: etaFD alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
					[netcdfFile waitUntilAllOperationsAreFinished];
					[rvHistory concatenateWithLowerDimensionalVariable: rv2 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
					[netcdfFile waitUntilAllOperationsAreFinished];
				}
				
				[xPositionHistory concatenateWithLowerDimensionalVariable: xPosition alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[netcdfFile waitUntilAllOperationsAreFinished];
				[yPositionHistory concatenateWithLowerDimensionalVariable: yPosition alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[netcdfFile waitUntilAllOperationsAreFinished];
			}
		}
		
		NSLog(@"Close the NetCDF file and wrap up");
		
		[equation waitUntilAllOperationsAreFinished];
		
		// The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
		[netcdfFile waitUntilAllOperationsAreFinished];
		[netcdfFile close];
		
		
	}
	return 0;
}

