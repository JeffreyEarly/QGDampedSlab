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
        
        GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: qg.dimensions forEquation: equation];
        
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
        
#warning For some reason I get total crap if I don't convert from a NetCDFVariable to a regular variable.
        u_wind = [GLFunction functionFromFunction: u_wind];
        v_wind = [GLFunction functionFromFunction: v_wind];
        
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
        /*		Create a differential operator for damping the slab layer                               */
        /************************************************************************************************/
        
        GLFloat k = pow(-1, dampingOrder+1)*pow(xDim.sampleInterval/M_PI,2*dampingOrder)/dampingTime;
        
        NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
        GLLinearTransform *harmonicOp = [GLLinearTransform harmonicOperatorOfOrder: dampingOrder fromDimensions: spectralDimensions forEquation: equation];
        GLLinearTransform *dampingOp = [[harmonicOp times: @(k)] plus: @(-1/linearDampingTime)];
		
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
        
        GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [[NSURL fileURLWithPath: [NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) firstObject]] URLByAppendingPathComponent:@"QGDampedSlab.nc"] forEquation: equation overwriteExisting: YES];
        
        GLFunction *dimensionalEta1 = [eta1_0 scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
        GLMutableVariable *eta1History = [dimensionalEta1 variableByAddingDimension: tDim];
        eta1History.name = @"eta1";
        eta1History.units = @"m";
        eta1History = [netcdfFile addVariable: eta1History];
        
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
		
        GLFunction *dimensionalXPosition = [xPosition scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *xPositionHistory = [dimensionalXPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
		
        GLFunction *dimensionalYPosition = [yPosition scaleVariableBy: qg.L_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
		GLMutableVariable *yPositionHistory = [dimensionalYPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
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
        GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[uI_0, vI_0, eta1_0_FD] stepSize: timeStep fFromTY:^(GLScalar *t, NSArray *yNew) {
            GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[t]];
            GLFunction *tx = interp.result[0];
            GLFunction *ty = interp.result[1];
            
            GLFunction *uI = yNew[0];
            GLFunction *vI = yNew[1];
            GLFunction *eta1 = yNew[2];
            
            // Should I be damping uI? Or u1?
            GLFunction *f_uI = [[[vI times:  @(u_scale)] plus: tx] plus: [[[[eta1 x] times: @(u_scale)] plus: [uI differentiateWithOperator:dampingOp]] spatialDomain]];
            GLFunction *f_vI = [[[[uI negate] times:  @(u_scale)] plus: ty] plus: [[[[eta1 y] times: @(u_scale)] plus: [vI differentiateWithOperator:dampingOp]] spatialDomain]];
            GLFunction *f_eta1 = [[[uI x] plus: [vI y]] times: @(div_scale)];
            
            NSArray *f = @[f_uI, f_vI, f_eta1];
            return f;
        }];
		
        
        for (GLFloat time = 0; time < maxTime/qg.T_QG; time += (86400/20)/qg.T_QG)
        {
            @autoreleasepool {
                NSArray *yin = [integrator stepForwardToTime: time];
                
                NSLog(@"Logging day: %f, step size: %f.", (qg.T_QG*integrator.currentTime/86400), integrator.lastStepSize*qg.T_QG);
                
                [tDim addPoint: @(time*qg.T_QG)];
                
                GLFunction *eta1 = [[yin[2] spatialDomain] scaleVariableBy: qg.N_QG withUnits: @"m" dimensionsBy: qg.L_QG units: @"m"];
                [eta1History concatenateWithLowerDimensionalVariable: eta1 alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLFunction *u = [yin[0] scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
                [uHistory concatenateWithLowerDimensionalVariable: u alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLFunction *v = [yin[1] scaleVariableBy: U_scale withUnits: @"m/s" dimensionsBy: qg.L_QG units: @"m"];
                [vHistory concatenateWithLowerDimensionalVariable: v alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[[GLScalar scalarWithValue: time forEquation:equation]]];
                
                GLScalar *tau_x_i = [interp.result[0] scaleBy: 1/tau_scale withUnits: @"N/m^2"];
                GLScalar *tau_y_i = [interp.result[1] scaleBy: 1/tau_scale withUnits: @"N/m^2"];
                [tau_x_i solve];
                [tau_y_i solve];
                [tau_xHistory concatenateWithLowerDimensionalVariable: tau_x_i alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [tau_yHistory concatenateWithLowerDimensionalVariable: tau_y_i alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                
                NSLog(@"Tau (%g, %g)", *(tau_x_i.pointerValue), *(tau_y_i.pointerValue));
                
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

