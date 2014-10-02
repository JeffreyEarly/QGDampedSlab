//
//  main.swift
//  DampedSlabSwift
//
//  Created by Jeffrey J. Early on 9/30/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

import Foundation
import GLNumericalModelingKit

let H0 = 100.0
let dRho1 = 1E-4
let maxTime = 40.0*86400.0
let dampingOrder=2.0 // order of the damping operator. Order 1 is harmonic, order 2 is biharmonic, etc.
let dampingTime=3600.0 // e-folding time scale of the Nyquist frequency.
let linearDampingTime = 4.0*86400.0

let f = 2 * 7.2921E-5 * sin( 24*M_PI/180 )
let g = 9.81
let g1 = dRho1 * g
let rho_water = 1025.0

let xDim = GLDimension(dimensionWithGrid: .PeriodicGrid, nPoints: 256, domainMin: -50e3, length: 100e3)
xDim.name = "x"
xDim.units = "m"

let yDim = GLDimension(dimensionWithGrid: .PeriodicGrid, nPoints: 256, domainMin: -50e3, length: 100e3)
yDim.name = "y"
yDim.units = "m"

let tDim = GLMutableDimension(points: [0.0])
tDim.name = "time"
tDim.units  = "s"

let equation = GLEquation()

let spatialDimensions = [xDim,yDim]
let x = GLFunction(ofRealTypeFromDimension: xDim, withDimensions: spatialDimensions, forEquation: equation)
let y = GLFunction(ofRealTypeFromDimension: yDim, withDimensions: spatialDimensions, forEquation: equation)

var h0 = GLFunction(ofRealTypeWithDimensions: spatialDimensions, forEquation: equation)
let u0 = GLFunction(ofRealTypeWithDimensions: spatialDimensions, forEquation: equation)
let v0 = GLFunction(ofRealTypeWithDimensions: spatialDimensions, forEquation: equation)

h0 = h0.setValue(H0, atIndices: ":,:");
u0.zero()
v0.zero()

println("Spatial resolution: \(xDim.sampleInterval) m")



