type Vector =  Array<number>
type Matrix = Array<Vector>
type Tensor = Array<Matrix>
type matrixRange = [number, number]


interface BoundaryCondition<T> {
  R1: T
  R2: T
  deltaTime: T
  rho: T
  Crho: T
  lambda: T
  range: T
  T0: Array<T>
  alpha1: T
  alpha2: T
  T1: T
  T2: T
}

interface Abbreviation<T> {
  a: T
  k1: T
  k2: T
  deltaT: T
  deltaX: T
}

let boundaryCondition: BoundaryCondition<number> = {
  R1: 10, //m
  R2: 2000, //m
  deltaTime: 10000, //s
  rho: 1000, //kg/m^3
  Crho: 100, //J/kg
  lambda: 100, //W/(m*K)
  range: 5, //T0.length
  T0: [50, 0, 0, 0, 300], 
  alpha1: 10, //W/(m^2*K)
  alpha2: 50, //W/(m^2*K)
  T1: 10000, //K
  T2: 0, //K
  }

//Variant 24
const testMatrix: Matrix = [
  [8, 7, -3, 2],
  [6, -13, -8, -5],
  [1, 2, 1, 26],
  [-10, -1, -2, 2]
]
const testRightPart: Matrix = toMatrix([47, -61, 95, -70])
const testAnswer: Matrix = [[6], [2], [7], [3]]

//Variant 25
const testMatrix1: Matrix = [
  [20, 0, 4, -5],
  [0, -20, 3, 1],
  [-2, 5, 60, 3],
  [7, -2, 4, -20]
]
const testRightPart1: Matrix = toMatrix([13, -12, -2, 3])
const testAnswer1: Matrix = toMatrix([1, 3, -1, -3])

function abbr(bounds: BoundaryCondition<number>): Abbreviation<number> {
  const a: number = bounds.lambda / (bounds.Crho * bounds.rho)
  const k1: number = bounds.alpha1 / (bounds.Crho * bounds.rho)
  const k2: number = bounds.alpha2 / (bounds.Crho * bounds.rho)
  const deltaT: number = bounds.deltaTime
  const deltaX: number = (bounds.R2 - bounds.R1) / bounds.T0.length

  return {a, k1, k2, deltaT, deltaX}
}

function toMatrix(a: Vector): Matrix {
  let result = [Array(a.length)]
  for (let i = 0; i < a.length; i++) {
    result[i] = Array(1)
    result[i][0]= a[i];
  }
  return result
}

function matrixMultiplication(a: Matrix, b: Matrix): Matrix {
  //define ranges of elements
  const aRange = [a.length, (a[1] as any)?.length || 1]    
  const bRange = [b.length, (b[1] as any)?.length || 1]
  if (aRange[1] != bRange[1]) {
    throw new Error("Invalid operation")
  }
  let result: Matrix = Array(a.length)
  //result.fill(0,a.length,)
  //calculation
  for (let i = 0; i < a.length; i++) {
    result[i] = Array(b.length)
    for (let j = 0; j < b.length; j++) {
      result[i][j] = 0
      for (let l = 0; l < aRange[1]; l++) {
        result[i][j] += a[i][l] * b[l][j]
      }
    }
  }
  //resulting
  return result
}

function zeidel(initialMatrix: Matrix, rightPart: Matrix, accuracy: number): Matrix {
  let tempX: number = 0
  let matrixNorm: number
  const matrixRange: number = initialMatrix.length
  const result: Matrix = toMatrix(Array(matrixRange))
  const maxIterations = 10
  
  //initial
  for (let i = 0; i < matrixRange; i++) {
    result[i][0] = (rightPart[i][0])/(initialMatrix[i][i])  
  }
  console.log("result:", result);
  for (let iteration = 0; iteration <= maxIterations; iteration++) {
    matrixNorm = 0
    for (let i = 0; i < matrixRange; i++) {
      tempX = 0
      for (let j = 0; j < matrixRange; j++) {
        tempX += initialMatrix[i][j]*result[j][0]
      }
    tempX = result[i][0] + (rightPart[i][0] - tempX)/initialMatrix[i][i]
    if (Math.abs(result[i][0]-tempX) > matrixNorm) {
      matrixNorm = Math.abs(result[i][0]-tempX)
    }
    result[i][0] = tempX
    }
    console.log(result);
    if (matrixNorm < accuracy) {
      return result
    }
  }

  throw new Error("Bad");
  
  
}

function implicitSchema(initialTemperature: Vector, boundaryCondition: BoundaryCondition<number> ): Vector {
  const bc = boundaryCondition
  console.log('bc: ', bc);
  const inT = initialTemperature
  console.log('inT: ', inT);
  const {a, k1, k2, deltaT, deltaX}: Abbreviation<number> = abbr(bc)
  const r = (i: number): number => bc.R1 + (i - 1) * deltaX
  //init resulted vector with unphysical
  const finalT: Vector = new Array<number>(inT.length + 2)
  //first value adding in array
  const T00 = inT[2] + 2 * k1 * deltaX * (bc.T1 - inT[1]) //(inT[0] + k1 * deltaX * bc.T1) / (1 + k1 * deltaX)
  inT.unshift(T00)
  //last value adding in array
  const T0n = inT[inT.length - 2] + 2 * k2 * deltaX * (bc.T2 - inT[inT.length - 1]) //(inT[inT.length - 1] + k2 * deltaX * bc.T2) / (1 + k2 * deltaX)
  inT.push(T0n)
  console.log('inT after adding: ', inT);
  //finding in range of [1, n-1]
  for (let i = 1; i < inT.length - 1; i++) {
    finalT[i] = inT[i] + a * deltaT / (deltaX ^ 2) * (inT[i + 1] - 2 * inT[i] + inT[i - 1]) + 
    a * deltaT / (2 * deltaX * r(i)) * (inT[i + 1] - inT[i - 1])
  }

  return finalT//.slice(1, -1)
}

//console.log(zeidel(testMatrix1, testRightPart1, 0.001))
function makeMatrixFromBC(timeSteps: number, initialConditions: BoundaryCondition<number> = boundaryCondition ): any {
  const spaceSteps: number = initialConditions.T0.length
  const equationNumbers: number = spaceSteps * timeSteps
  // const coloumn: number[] = new Array(spaceSteps)
  // const matrix: Matrix = new Array<>
  let result: Array<Array<Array<number>>> = Array<Matrix>(equationNumbers)
  for (let equationNumber = 0; equationNumber < equationNumbers; equationNumber++) {
    result[equationNumber] = Array<Array<number>>(spaceSteps)
    for (let spaceStep = 0; spaceStep < spaceSteps; spaceStep++) {
      result[equationNumber][spaceStep] = Array<number>(timeSteps)
      for (let timeStep = 0; timeStep < timeSteps; timeStep++) {
        result[equationNumber][spaceStep][timeStep] = 0
      }
    }
  }
  console.log(result);
  return result
}


//исполняющий код
// let temp = makeMatrixFromBC(2, boundaryCondition)
// temp[0][1][1] = 1
// console.log('temp[0][1][2]: ', temp[0][1][1]);
// console.log(temp);
console.log(implicitSchema(boundaryCondition.T0, boundaryCondition))
