"use strict";
let boundaryCondition = {
    R1: 10,
    R2: 2000,
    deltaTime: 100,
    rho: 1000,
    Crho: 100,
    lambda: 100,
    range: 5,
    T0: [50, 100, 150, 200, 300],
    alpha1: 10,
    alpha2: 50,
    T1: 50,
    T2: 1000, //K
};
//Variant 24
const testMatrix = [
    [8, 7, -3, 2],
    [6, -13, -8, -5],
    [1, 2, 1, 26],
    [-10, -1, -2, 2]
];
const testRightPart = toMatrix([47, -61, 95, -70]);
const testAnswer = [[6], [2], [7], [3]];
//Variant 25
const testMatrix1 = [
    [20, 0, 4, -5],
    [0, -20, 3, 1],
    [-2, 5, 60, 3],
    [7, -2, 4, -20]
];
const testRightPart1 = toMatrix([13, -12, -2, 3]);
const testAnswer1 = toMatrix([1, 3, -1, -3]);
function abbr(bounds) {
    const a = bounds.lambda / (bounds.Crho * bounds.rho);
    const k1 = bounds.alpha1 / (bounds.Crho * bounds.rho);
    const k2 = bounds.alpha2 / (bounds.Crho * bounds.rho);
    const deltaT = bounds.deltaTime;
    const deltaX = (bounds.R2 - bounds.R1) / bounds.T0.length;
    return { a, k1, k2, deltaT, deltaX };
}
function toMatrix(a) {
    let result = [Array(a.length)];
    for (let i = 0; i < a.length; i++) {
        result[i] = Array(1);
        result[i][0] = a[i];
    }
    return result;
}
function matrixMultiplication(a, b) {
    var _a, _b;
    //define ranges of elements
    const aRange = [a.length, ((_a = a[1]) === null || _a === void 0 ? void 0 : _a.length) || 1];
    const bRange = [b.length, ((_b = b[1]) === null || _b === void 0 ? void 0 : _b.length) || 1];
    if (aRange[1] != bRange[1]) {
        throw new Error("Invalid operation");
    }
    let result = Array(a.length);
    //result.fill(0,a.length,)
    //calculation
    for (let i = 0; i < a.length; i++) {
        result[i] = Array(b.length);
        for (let j = 0; j < b.length; j++) {
            result[i][j] = 0;
            for (let l = 0; l < aRange[1]; l++) {
                result[i][j] += a[i][l] * b[l][j];
            }
        }
    }
    //resulting
    return result;
}
function zeidel(initialMatrix, rightPart, accuracy) {
    let tempX = 0;
    let matrixNorm;
    const matrixRange = initialMatrix.length;
    const result = toMatrix(Array(matrixRange));
    const maxIterations = 10;
    //initial
    for (let i = 0; i < matrixRange; i++) {
        result[i][0] = (rightPart[i][0]) / (initialMatrix[i][i]);
    }
    console.log("result:", result);
    for (let iteration = 0; iteration <= maxIterations; iteration++) {
        matrixNorm = 0;
        for (let i = 0; i < matrixRange; i++) {
            tempX = 0;
            for (let j = 0; j < matrixRange; j++) {
                tempX += initialMatrix[i][j] * result[j][0];
            }
            tempX = result[i][0] + (rightPart[i][0] - tempX) / initialMatrix[i][i];
            if (Math.abs(result[i][0] - tempX) > matrixNorm) {
                matrixNorm = Math.abs(result[i][0] - tempX);
            }
            result[i][0] = tempX;
        }
        console.log(result);
        if (matrixNorm < accuracy) {
            return result;
        }
    }
    throw new Error("Bad");
}
function implicitSchema(initialTemperature, boundaryCondition) {
    const bc = boundaryCondition;
    console.log('bc: ', bc);
    const inT = initialTemperature;
    console.log('inT: ', inT);
    const { a, k1, k2, deltaT, deltaX } = abbr(bc);
    const r = (i) => bc.R1 + (i - 1) * deltaX;
    //init resulted vector with unphysical
    const finalT = new Array(inT.length + 2);
    //first value adding in array
    const T00 = (inT[0] + k1 * deltaX * bc.T1) / (1 + k1 * deltaX);
    inT.unshift(T00);
    //last value adding in array
    const T0n = (inT[inT.length - 1] + k2 * deltaX * bc.T2) / (1 + k2 * deltaX);
    inT.push(T0n);
    console.log('inT after adding: ', inT);
    //finding in range of [1, n-1]
    for (let i = 1; i < inT.length - 1; i++) {
        finalT[i] = inT[i] + a * deltaT / (deltaX ^ 2) * (inT[i + 1] - 2 * inT[i] + inT[i - 1]) +
            a * deltaT / (2 * deltaX * r(i)) * (inT[i + 1] - inT[i - 1]);
    }
    return finalT; //.slice(1, -1)
}
//console.log(zeidel(testMatrix1, testRightPart1, 0.001))
function makeMatrixFromBC(timeSteps, initialConditions = boundaryCondition) {
    const spaceSteps = initialConditions.T0.length;
    const equationNumbers = spaceSteps * timeSteps;
    // const coloumn: number[] = new Array(spaceSteps)
    // const matrix: Matrix = new Array<>
    let result = Array(equationNumbers);
    for (let equationNumber = 0; equationNumber < equationNumbers; equationNumber++) {
        result[equationNumber] = Array(spaceSteps);
        for (let spaceStep = 0; spaceStep < spaceSteps; spaceStep++) {
            result[equationNumber][spaceStep] = Array(timeSteps);
            for (let timeStep = 0; timeStep < timeSteps; timeStep++) {
                result[equationNumber][spaceStep][timeStep] = 0;
            }
        }
    }
    console.log(result);
    return result;
}
//исполняющий код
// let temp = makeMatrixFromBC(2, boundaryCondition)
// temp[0][1][1] = 1
// console.log('temp[0][1][2]: ', temp[0][1][1]);
// console.log(temp);
console.log(implicitSchema(boundaryCondition.T0, boundaryCondition));
