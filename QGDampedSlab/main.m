//
//  main.m
//  QGDampedSlab
//
//  Created by Jeffrey J. Early on 10/17/12.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main (int argc, const char * argv[])
{
	
	@autoreleasepool {
		
		GLFloat H1 = 50;
		GLFloat H2 = 800;
		GLFloat dRho1 = 1E-3;
		GLFloat dRho2 = 1E-3;
		GLFloat latitude = 24;
		
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
		GLFloat L_1 = sqrt( gprime * H1) / f0;
		GLFloat L_2 = sqrt( gprime * H2) / f0;
		GLFloat u_scale = f0/(beta*L_2);
		GLFloat tau_scale = 1./( rho_water*H1*beta*beta*L_2*L_2*L_2);
		GLFloat div_scale = u_scale*(L_1*L_1)/(L_2*L_2);
		
		GLEquation *equation = [[GLEquation alloc] init];
		GLNetCDFFile *restartFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Volumes/Data/ForcedDissipativeQGTurbulence/QGTurbulenceRestart_kFrac_16.nc"] forEquation: equation];
		
		GLFloat k_f = [restartFile.globalAttributes[@"forcing_wavenumber"] doubleValue];
		GLFloat forcingWidth = [restartFile.globalAttributes[@"forcing_width"] doubleValue];
		GLFloat nu = [restartFile.globalAttributes[@"nu"] doubleValue];
		GLFloat alpha = [restartFile.globalAttributes[@"alpha"] doubleValue];
		GLFloat energy = [restartFile.globalAttributes[@"energy"] doubleValue];
		BOOL shouldUseSVV = [restartFile.globalAttributes[@"uses-spectral-vanishing-viscosity"] boolValue];
		
		[GLBasisTransformOperation setShouldAntialias: [restartFile.globalAttributes[@"is-anti-aliased"] boolValue]];
		
		GLFloat T_QG = [restartFile.globalAttributes[@"time_scale"] doubleValue] /86400; // days
		
		maxTime = 10;
		
		BOOL writeExtrasToFile = NO;
		
		/************************************************************************************************/
		/*		Read the winds from file																*/
		/************************************************************************************************/
		
		// If all goes well, the variable t will be identified as the coordinated variable and therefore turned into a dimension, leaving only u and v.
		GLNetCDFFile *winds = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Volumes/Data/QGPlusSlab/winds.nc"] forEquation: equation];
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
		
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = restartFile.dimensions[0];
		GLDimension *yDim = restartFile.dimensions[1];
		
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		[equation setDefaultDifferentiationBasis: @[@(kGLExponentialBasis), @(kGLExponentialBasis)] forOrder: 2];
		
		GLVariable *ssh = restartFile.variables[0];
		[ssh solve];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLVariable *x = [GLVariable variableOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLVariable *y = [GLVariable variableOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		// At the moment we know that this is the spectral operators, although in the future we'll have to set this up explicitly.
		GLSpectralDifferentialOperatorPool *diffOperators = [equation defaultDifferentialOperatorPoolForVariable: ssh];
		
		// Create the operator xx+yy-1---this is how you compute y from eta
		GLSpectralDifferentialOperator *laplacianMinusOne = [[diffOperators harmonicOperator] scalarAdd: -1.0];
		[diffOperators setDifferentialOperator: laplacianMinusOne forName: @"laplacianMinusOne"];
		
		// Create the operator 1/(xx+yy-1)---this is how you compute eta from y.
		[diffOperators setDifferentialOperator: [laplacianMinusOne scalarDivide: 1.0] forName: @"inverseLaplacianMinusOne"];
		
		// This builds the differentiation matrix diff_{xxx} + diff_{xyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxx] plus: [diffOperators xyy]] forName: @"diffJacobianX"];
		
		// This builds the differentiation matrix diff_{xxy} + diff_{yyy}
		[diffOperators setDifferentialOperator: [[diffOperators xxy] plus: [diffOperators yyy]] forName: @"diffJacobianY"];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		GLSpectralDifferentialOperator *viscosity;
		if (shouldUseSVV) {
			GLSpectralDifferentialOperator *svv = [diffOperators spectralVanishingViscosityFilter];
			viscosity = [[[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: nu] multiply: svv];
		} else {
			viscosity = [[diffOperators harmonicOperatorOfOrder: 2] scalarMultiply: nu];
		}
		
		viscosity = [viscosity scalarAdd: alpha];
		
		[diffOperators setDifferentialOperator: viscosity forName: @"diffLin"];
		
		
		/************************************************************************************************/
		/*		Plop down a float at each grid point													*/
		/************************************************************************************************/
		
		GLVariable *xPosition = [GLVariable variableFromVariable: x];
		GLVariable *yPosition = [GLVariable variableFromVariable: y];
		
		GLVariable *sshFD = [ssh frequencyDomain];
		GLVariable *rv = [[ssh diff:@"harmonicOperator"] spatialDomain];
		
		/************************************************************************************************/
		/*		Create a NetCDF file to record everything interesting.									*/
		/************************************************************************************************/
		
		// Now we create a mutable variable in order to record the evolution of the Gaussian.
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: @"/Volumes/Data/QGTurbulenceRestarted_2.nc"] forEquation: equation overwriteExisting: YES];
		for (id key in restartFile.globalAttributes) {
			[netcdfFile setGlobalAttribute: restartFile.globalAttributes[key] forKey: key];
		}
		
		GLMutableVariable *sshHistory = [ssh variableByAddingDimension: tDim];
		sshHistory.name = @"SSH";
		sshHistory = [netcdfFile addVariable: sshHistory];
		
		GLMutableVariable *sshFDHistory, *rvHistory;
		if (writeExtrasToFile) {
			sshFDHistory = [sshFD variableByAddingDimension: tDim];
			sshFDHistory.name = @"SSH_FD";
			sshFDHistory = [netcdfFile addVariable: sshFDHistory];
			
			rvHistory = [rv variableByAddingDimension: tDim];
			rvHistory.name = @"RV";
			rvHistory = [netcdfFile addVariable: rvHistory];
		}
		
		GLMutableVariable *xPositionHistory = [xPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
		
		GLMutableVariable *yPositionHistory = [yPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
		
		/************************************************************************************************/
		/*		Create an integrator for the qg equation												*/
		/************************************************************************************************/
		
		CGFloat cfl = 0.5;
		GLFloat viscous_dt = 0.1 * xDim.sampleInterval * xDim.sampleInterval / nu;
		
		GLVariable *zeta = [ssh diff: @"laplacianMinusOne"];
		NSArray *yin = @[zeta, xPosition, yPosition];
		
		GLAdaptiveRungeKuttaOperation *qgIntegrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: yin stepSize: viscous_dt fFromTY:^(GLVariable *time, NSArray *yNew) {
			GLVariable *eta = [yNew[0] diff: @"inverseLaplacianMinusOne"];
			GLVariable *f = [[eta diff:@"diffLin"] plus: [[[[eta y] times: [eta diff: @"diffJacobianX"]] minus: [[eta x] times: [eta diff: @"diffJacobianY"]]] frequencyDomain]];
			
			NSArray *uv = @[[[[eta y] spatialDomain] negate], [[eta x] spatialDomain] ];
			NSArray *xy = @[yNew[1], yNew[2]];
			GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
			
			return @[f, interp.result[0], interp.result[1]];
		}];
		
		/************************************************************************************************/
		/*		Determine an appropriate time step based on the CFL condition.							*/
		/************************************************************************************************/
		
		GLFloat inertialTimeStep = 600/T_QG;
		
		GLVectorIntegrationOperation *integrator = [GLVectorIntegrationOperation rungeKutta4AdvanceY: yin stepSize: inertialTimeStep fFromTY:^(GLVariable *t, NSArray *yNew) {
			
			GLInterpolationOperation *interp = [[GLInterpolationOperation alloc] initWithFirstOperand: @[tau_x, tau_y] secondOperand: @[[t scalarMultiply: T_QG*86400]]];
			GLVariable *tx = interp.result[0];
			GLVariable *ty = interp.result[1];
			
			GLVariable *u = yNew[0];
			GLVariable *v = yNew[1];
			GLVariable *h = yNew[2];
			
			//			GLVariable *fu = [[[[v scalarMultiply: f] plus: [[tx scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h x] scalarMultiply: g1] minus: [u diff:@"damp"]] spatialDomain]] minus: [[u times: [u x]] plus: [v times: [u y]]]];
			//            GLVariable *fv = [[[[u scalarMultiply: -f] plus: [[ty scalarMultiply: 1/rho_water] dividedBy: h]] minus: [[[[h y] scalarMultiply: g1] minus: [v diff:@"damp"]] spatialDomain]] minus: [[u times: [v x]] plus: [v times: [v y]]]];
			
			GLVariable *fu = [[[v scalarMultiply: u_scale] plus: [tx scalarMultiply: tau_scale]] minus: [[[[h x] scalarMultiply: u_scale] minus: [u diff:@"damp"]] spatialDomain]];
			GLVariable *fv = [[[u scalarMultiply: -u_scale] plus: [ty scalarMultiply: tau_scale]] minus: [[[[h y] scalarMultiply: u_scale] minus: [v diff:@"damp"]] spatialDomain]];
			GLVariable *fh = [[[[u x] plus: [v y]] scalarMultiply:div_scale] negate];
			
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

