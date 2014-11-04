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
        NSUInteger floatFraction = 4;
        
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
        GLFloat U_scale = beta*L_2*L_2;
		
        
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
                
        GLFunction *eta1_0 = [GLFunction functionOfRealTypeWithDimensions: qg.dimensions forEquation: equation];
        GLFunction *uI_0 = [GLFunction functionOfRealTypeWithDimensions: qg.dimensions forEquation: equation];
        GLFunction *vI_0 = [GLFunction functionOfRealTypeWithDimensions: qg.dimensions forEquation: equation];
        
        [eta1_0 zero];
        [uI_0 zero];
        [vI_0 zero];
		
        /************************************************************************************************/
        /*		Read the winds from file																*/
        /************************************************************************************************/
        
        // If all goes well, the variable t will be identified as the coordinated variable and therefore turned into a dimension, leaving only u and v.
        GLNetCDFFile *winds = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Users/jearly/Documents/Models/QGDampedSlab/winds.nc"] forEquation: equation];
        GLFunction *u_wind = winds.variables[0];
        GLFunction *v_wind = winds.variables[1];
		
        GLFloat rho_air = 1.25;
        GLFunction *speed_wind = [[[u_wind times: u_wind] plus: [v_wind times: v_wind]] sqrt];
        GLFunction *dragCoefficient = [[[speed_wind scalarMultiply: 0.071] scalarAdd: 0.50] scalarMultiply: 1e-3];
        GLFunction *tau_x = [[[u_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air];
        GLFunction *tau_y = [[[v_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air];
        
        GLFloat tau_scale = 1./( rho_water*H1*beta*beta*L_2*L_2*L_2); // Nondimensionalization factor
        tau_x = [tau_x scaleVariableBy: tau_scale withUnits: @"unitless" dimensionsBy:1/qg.T_QG units: @"unitless"];
        tau_y = [tau_y scaleVariableBy: tau_scale withUnits: @"unitless" dimensionsBy:1/qg.T_QG units: @"unitless"];
        
        [tau_x solve];
        [tau_y solve];
		
		/************************************************************************************************/
		/*		Plop down floats                                                                        */
		/************************************************************************************************/
		
		GLDimension *xDimND = qg.dimensions[0];
		GLDimension *yDimND = qg.dimensions[1];
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: xDimND.gridType nPoints: xDimND.nPoints/floatFraction domainMin: xDimND.domainMin length:xDimND.domainLength];
        xFloatDim.name = @"x-float";
        GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: yDimND.gridType nPoints: yDimND.nPoints/floatFraction domainMin: yDimND.domainMin length:yDimND.domainLength];
        yFloatDim.name = @"y-float";
        
        NSArray *floatDims = @[xFloatDim, yFloatDim];
        GLFunction *xPos1 = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDims forEquation: equation];
        GLFunction *yPos1 = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDims forEquation: equation];
		GLFunction *xPos2 = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDims forEquation: equation];
		GLFunction *yPos2 = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDims forEquation: equation];
		
        /************************************************************************************************/
        /*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
        /************************************************************************************************/
        
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [[NSURL fileURLWithPath: [NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) firstObject]] URLByAppendingPathComponent:@"QGDampedSlab.nc"] forEquation: equation overwriteExisting: YES];
        
        GLFunction *dimensionalEta1 = [eta1_0 scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
        GLMutableVariable *eta1History = [dimensionalEta1 variableByAddingDimension: tDim];
        eta1History.name = @"eta1";
        eta1History.units = @"m";
        eta1History = [netcdfFile addVariable: eta1History];
        
        GLFunction *dimensionalEta2 = [[qg.ssh spatialDomain] scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
        GLMutableVariable *eta2History = [dimensionalEta2 variableByAddingDimension: tDim];
        eta2History.name = @"eta2";
        eta2History.units = @"m";
        eta2History = [netcdfFile addVariable: eta2History];
        
        GLFunction *dimensionalU = [uI_0 scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
        GLMutableVariable *uHistory = [dimensionalU variableByAddingDimension: tDim];
        uHistory.name = @"u";
        uHistory.units = @"m/s";
        uHistory = [netcdfFile addVariable: uHistory];
        
        GLFunction *dimensionalV = [vI_0 scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
        GLMutableVariable *vHistory = [dimensionalV variableByAddingDimension: tDim];
        vHistory.name = @"v";
        vHistory.units = @"m/s";
        vHistory = [netcdfFile addVariable: vHistory];
        
        GLFunction *tau0_x = [GLFunction functionOfRealTypeWithDimensions: @[] forEquation: equation];
		GLMutableVariable *tau_xHistory = [tau0_x variableByAddingDimension: tDim];
		tau_xHistory.name = @"tau_x";
        tau_xHistory.units = @"N/m^2";
		tau_xHistory = [netcdfFile addVariable: tau_xHistory];
        
        GLFunction *tau0_y = [GLFunction functionOfRealTypeWithDimensions: @[] forEquation: equation];
		GLMutableVariable *tau_yHistory = [tau0_y variableByAddingDimension: tDim];
		tau_yHistory.name = @"tau_y";
        tau_yHistory.units = @"N/m^2";
		tau_yHistory = [netcdfFile addVariable: tau_yHistory];
		
        GLFunction *dimensionalXPosition1 = [xPos1 scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *xPosition1History = [dimensionalXPosition1 variableByAddingDimension: tDim];
		xPosition1History.name = @"x-position-layer-1";
		xPosition1History = [netcdfFile addVariable: xPosition1History];
		
        GLFunction *dimensionalYPosition1 = [yPos1 scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *yPosition1History = [dimensionalYPosition1 variableByAddingDimension: tDim];
		yPosition1History.name = @"y-position-layer-1";
		yPosition1History = [netcdfFile addVariable: yPosition1History];
		
		GLFunction *dimensionalXPosition2 = [xPos2 scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *xPosition2History = [dimensionalXPosition2 variableByAddingDimension: tDim];
		xPosition2History.name = @"x-position-layer-2";
		xPosition2History = [netcdfFile addVariable: xPosition2History];
		
		GLFunction *dimensionalYPosition2 = [yPos2 scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *yPosition2History = [dimensionalYPosition2 variableByAddingDimension: tDim];
		yPosition2History.name = @"y-position-layer-2";
		yPosition2History = [netcdfFile addVariable: yPosition2History];
		
        /************************************************************************************************/
        /*		Determine an appropriate time step based on the CFL condition.							*/
        /************************************************************************************************/
        
        CGFloat cfl = 0.5;
        GLFloat timeStep = cfl * xDim.sampleInterval / sqrt(g*H1);
        timeStep = 300/qg.T_QG;
		
		/************************************************************************************************/
		/*		Determine an appropriate time step based on the CFL condition.							*/
		/************************************************************************************************/
        
        GLFloat u_scale = f0/(beta*L_2);
        
        GLFloat div_scale = u_scale*(L_1*L_1)/(L_2*L_2);
        GLFunction *eta1_0_FD = [eta1_0 transformToBasis: @[@(kGLExponentialBasis),@(kGLExponentialBasis)]];
        [qg createDifferentialOperators];
        [qg createIntegrationOperation];
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[qg.yin[0], uI_0, vI_0, eta1_0_FD, xPos1, yPos1, xPos2, yPos2] stepSize: timeStep fFromTY:^(GLScalar *t, NSArray *y) {
            // This should work because the first item in y is expected to be the same for both blocks (this one and the qg one).
            NSArray *fqg = qg.fBlock(t,y);
            GLFunction *eta2 = y[0];
            
            GLSimpleInterpolationOperation *interpStress = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[t]];
            GLFunction *tx = interpStress.result[0];
            GLFunction *ty = interpStress.result[1];
						
            GLFunction *uI = y[1];
            GLFunction *vI = y[2];
            GLFunction *eta1 = y[3];
            
            GLFunction *u2 = [[[eta2 y] spatialDomain] negate];
            GLFunction *v2 = [[eta2 x] spatialDomain];
            
            // Should I be damping uI? Or u1?
			// f_u = v + tx + eta_x + damp_u
			// f_v = -u + ty + eta_y + damp_v
			GLFunction *f_uI = [[[vI times:  @(u_scale)] plus: [tx minus: [uI times: @(qg.T_QG/linearDampingTime)]]] plus: [[[eta1 x] times: @(u_scale)] spatialDomain]];
            GLFunction *f_vI = [[[[uI negate] times:  @(u_scale)] plus: [ty minus: [vI times: @(qg.T_QG/linearDampingTime)]]] plus: [[[eta1 y] times: @(u_scale)] spatialDomain]];
            GLFunction *f_eta1 = [[[uI x] plus: [vI y]] times: @(div_scale)];
			
			GLSimpleInterpolationOperation *interpUV1 = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[uI,vI] secondOperand: @[y[4],y[5]]];
			
			GLSimpleInterpolationOperation *interpUV2 = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2,v2] secondOperand: @[y[5],y[6]]];
            
            NSArray *f = @[fqg[0], f_uI, f_vI, f_eta1, interpUV1.result[0], interpUV1.result[1], interpUV2.result[0], interpUV2.result[1]];
            return f;
        }];
		
        
        for (GLFloat time = (86400/20)/qg.T_QG; time < maxTime/qg.T_QG; time += (86400/20)/qg.T_QG)
        {
            @autoreleasepool {
                NSArray *yin = [integrator stepForwardToTime: time];
                
                NSLog(@"Logging day: %f, step size: %f.", (qg.T_QG*integrator.currentTime/86400), integrator.lastStepSize*qg.T_QG);
                
                [tDim addPoint: @(time*qg.T_QG)];
                
                GLFunction *eta2 = [[yin[0] spatialDomain] scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
                [eta2History concatenateWithLowerDimensionalVariable: eta2 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLFunction *u = [yin[1] scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
                [uHistory concatenateWithLowerDimensionalVariable: u alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLFunction *v = [yin[2] scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
                [vHistory concatenateWithLowerDimensionalVariable: v alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLFunction *eta1 = [[yin[3] spatialDomain] scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
                [eta1History concatenateWithLowerDimensionalVariable: eta1 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
				GLFunction *xPos1 = [yin[4] scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
				[xPosition1History concatenateWithLowerDimensionalVariable: xPos1 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
				GLFunction *yPos1 = [yin[5] scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
				[yPosition1History concatenateWithLowerDimensionalVariable: yPos1 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
				GLFunction *xPos2 = [yin[6] scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
				[xPosition2History concatenateWithLowerDimensionalVariable: xPos2 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
				GLFunction *yPos2 = [yin[7] scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
				[yPosition2History concatenateWithLowerDimensionalVariable: yPos2 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				
                GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[[GLScalar scalarWithValue: time forEquation:equation]]];
                
                GLScalar *tau_x_i = [interp.result[0] scaleBy: 1/tau_scale withUnits: @"N/m^2"];
                GLScalar *tau_y_i = [interp.result[1] scaleBy: 1/tau_scale withUnits: @"N/m^2"];
                [tau_x_i solve];
                [tau_y_i solve];
                [tau_xHistory concatenateWithLowerDimensionalVariable: tau_x_i alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [tau_yHistory concatenateWithLowerDimensionalVariable: tau_y_i alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
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

