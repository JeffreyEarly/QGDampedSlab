//
//  main.m
//  ShallowWaterWithWind
//
//  Created by Jeffrey J. Early on 10/11/12.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main(int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		GLFloat H0 = 100;
		GLFloat dRho1 = 1E-4;
		GLFloat maxTime = 40*86400;
		GLFloat dampingOrder=2; // order of the damping operator. Order 1 is harmonic, order 2 is biharmonic, etc.
		GLFloat dampingTime=3600; // e-folding time scale of the Nyquist frequency.
		GLFloat linearDampingTime = 4*86400;
		
		GLFloat f = 2 * 7.2921E-5 * sin( 24*M_PI/180. );
		GLFloat g = 9.81;
		GLFloat g1 = dRho1 * g;
		GLFloat rho_water = 1025;
		
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 256 domainMin: -50*1000 length: 100*1000];
		xDim.name = @"x";
		xDim.units = @"m";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: 256 domainMin: -50*1000 length: 100*1000];
		yDim.name = @"y";
		yDim.units = @"m";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		tDim.units = @"s";
		
		NSLog(@"Spatial resolution: %f m", xDim.sampleInterval);
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create the initial conditions															*/
		/************************************************************************************************/
		
		GLFunction *h0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
		GLFunction *u0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
		GLFunction *v0 = [GLFunction functionOfRealTypeWithDimensions: spatialDimensions forEquation: equation];
		
		h0 = [h0 setValue: H0 atIndices: @":,:"];
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
		GLFunction *speed_wind = [[[u_wind times: u_wind] plus: [v_wind times: v_wind]] sqrt];
		GLFunction *dragCoefficient = [[[speed_wind scalarMultiply: 0.071] scalarAdd: 0.50] scalarMultiply: 1e-3];
		GLFunction *tau_x = [[[u_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air];
		GLFunction *tau_y = [[[v_wind times: speed_wind] times: dragCoefficient] scalarMultiply: rho_air];
		
		[tau_x solve];
		[tau_y solve];
		
		//		GLVariable *tau0_x = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation:equation];
		//		GLVariable *tau0_y = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation:equation];
		//
		//		*(tau0_x.pointerValue) = *(tau_x.pointerValue);
		//		*(tau0_y.pointerValue) = *(tau_y.pointerValue);
		
		/************************************************************************************************/
		/*		Create and cache the differential operators that we will be using	                 	*/
		/************************************************************************************************/
		
		GLFloat k = pow(-1, dampingOrder+1)*pow(xDim.sampleInterval/M_PI,2*dampingOrder)/dampingTime;
		
		NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
		GLLinearTransform *harmonicOp = [GLLinearTransform harmonicOperatorOfOrder: dampingOrder fromDimensions: spectralDimensions forEquation: equation];
		GLLinearTransform *dampingOp = [[harmonicOp times: @(k)] plus: @(-1/linearDampingTime)];
		
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
		
		//		GLMutableVariable *tau_xHistory = [tau0_x variableByAddingDimension: tDim];
		//		tau_xHistory.name = @"tau_x";
		//        tau_xHistory.units = @"N/m^2";
		//		tau_xHistory = [netcdfFile addVariable: tau_xHistory];
		//
		//		GLMutableVariable *tau_yHistory = [tau0_y variableByAddingDimension: tDim];
		//		tau_yHistory.name = @"tau_y";
		//        tau_yHistory.units = @"N/m^2";
		//		tau_yHistory = [netcdfFile addVariable: tau_yHistory];
		
		/************************************************************************************************/
		/*		Determine an appropriate time step based on the CFL condition.							*/
		/************************************************************************************************/
		
		CGFloat cfl = 0.5;
		GLFloat timeStep = cfl * xDim.sampleInterval / sqrt(g1*H0);
		timeStep = 600;
		
		/************************************************************************************************/
		/*		Create the integration object.															*/
		/************************************************************************************************/
		
		//		GLVariable *tx = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: equation];
		//		GLVariable *ty = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: equation];
		//		*(tx.pointerValue) = 0.0;
		//		*(ty.pointerValue) = 0.5;
		
		
		//		GLVariable *u = u0;
		//		GLVariable *v = v0;
		//		GLVariable *h = h0;
		//
		//		GLVariable *fu = [[[v scalarMultiply: f] plus: [[tx scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h x] scalarMultiply: g1] minus: [u diff:@"damp"]] spatialDomain]];
		//		GLVariable *fv = [[[u scalarMultiply: -f] plus: [[ty scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h y] scalarMultiply: g1] minus: [v diff:@"damp"]] spatialDomain]];
		//		GLVariable *fh = [[[[u x] plus: [v y]] times: h] negate];
		//
		//		[fu solve]; [fv solve]; [fh solve];
		
		GLRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: @[u0, v0, h0] stepSize: timeStep fFromTY:^(GLScalar *t, NSArray *yNew) {
			
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand:  @[tau_x, tau_y] secondOperand: @[t]];
			GLFunction *tx = interp.result[0];
			GLFunction *ty = interp.result[1];
			
			GLFunction *u = yNew[0];
			GLFunction *v = yNew[1];
			GLFunction *h = yNew[2];
			
			//			GLVariable *fu = [[[[v scalarMultiply: f] plus: [[tx scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h x] scalarMultiply: g1] minus: [u diff:@"damp"]] spatialDomain]] minus: [[u times: [u x]] plus: [v times: [u y]]]];
			//            GLVariable *fv = [[[[u scalarMultiply: -f] plus: [[ty scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h y] scalarMultiply: g1] minus: [v diff:@"damp"]] spatialDomain]] minus: [[u times: [v x]] plus: [v times: [v y]]]];
			
			GLFunction *fu = [[[v times: @(f)] plus: [[tx times: @(1/rho_water)] dividedBy: h]] minus: [[[[h x] scalarMultiply: g1] minus: [u differentiateWithOperator:dampingOp]] spatialDomain]];
			GLFunction *fv = [[[u times: @(-f)] plus: [[ty times: @(1/rho_water)] dividedBy: h]] minus: [[[[h y] scalarMultiply: g1] minus: [v differentiateWithOperator:dampingOp]] spatialDomain]];
			GLFunction *fh = [[[[u x] plus: [v y]] times: h] negate];
			
			NSArray *f = @[fu, fv, fh];
			return f;
		}];
		
		/************************************************************************************************/
		/*		Step forward in time, and write the data to file every-so-often.						*/
		/************************************************************************************************/
		
		for (GLFloat time = 0; time < maxTime; time += 86400/20)
		{
			@autoreleasepool {
				NSArray *yin = [integrator stepForwardToTime: time];
				
				NSLog(@"Logging day: %f, step size: %f.", (integrator.currentTime/86400), integrator.lastStepSize);
				
				[tDim addPoint: @(time/86400)];
				
				[uHistory concatenateWithLowerDimensionalVariable: yin[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[vHistory concatenateWithLowerDimensionalVariable: yin[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
				[hHistory concatenateWithLowerDimensionalVariable: yin[2] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
			}
		}
		
		[equation waitUntilAllOperationsAreFinished];
		
		// The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
		[netcdfFile waitUntilAllOperationsAreFinished];
		[netcdfFile close];
	}
	return 0;
}
