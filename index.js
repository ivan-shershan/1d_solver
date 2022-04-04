"use strict";
let boundaryCondition = {
    R1: 10,
    R2: 2000,
    deltaTime: 10,
    rho: 1000,
    Crho: 100,
    lambda: 100,
    range: 5,
    T0: [50, 0, 0, 0, 300],
    alpha1: 10,
    alpha2: 50,
    T1: 10000,
    T2: 0, //K
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
function swipeMethod(matrix, rightPart) {
    const A = matrix;
    const b = rightPart;
    const n = b.length;
    //swipe coefficients
    let alpha = new Array(n);
    let beta = new Array(n);
    //first row
    let y = A[0][0];
    alpha[0] = -A[0][1] / y;
    beta[0] = b[0] / y;
    //other except last
    for (let i = 1; i < n - 1; i++) {
        y = A[i][i] + A[i][i - 1] * alpha[i - 1];
        alpha[i] = -A[i][i + 1] / y;
        beta[i] = (b[i] - A[i][i - 1] * beta[i - 1]) / y;
    }
    //last row
    y = A[n - 1][n - 1] + A[n - 1][n - 2] * alpha[n - 2];
    console.log("y", y);
    beta[n - 1] = (b[n - 1] - A[n - 1][n - 2] * beta[n - 2]) / y;
    console.log("alpha", alpha);
    console.log("beta", beta);
    //at last
    let answer = new Array(n);
    answer[n - 1] = beta[n - 1];
    for (let i = n - 2; i >= 0; i--) {
        answer[i] = alpha[i] * answer[i + 1] + beta[i];
    }
    return answer;
}
function swipeMethod1(a, b, c, d) {
    const n = b.length;
    //swipe coefficients
    let alpha = new Array(n);
    let beta = new Array(n);
    //first row
    let y = b[0];
    alpha[0] = -c[0] / y;
    beta[0] = d[0] / y;
    //other except last
    for (let i = 1; i < n - 1; i++) {
        y = b[i] + a[i] * alpha[i - 1];
        alpha[i] = -c[i] / y;
        beta[i] = (d[i] - a[i] * beta[i - 1]) / y;
    }
    //last row
    y = b[n - 1] + a[n - 1] * alpha[n - 2];
    console.log("y", y);
    beta[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 2]) / y;
    // console.log("abc", a,b,c)
    // console.log("alpha", alpha)
    // console.log("beta", beta)
    //at last
    let answer = new Array(n);
    answer[n - 1] = beta[n - 1];
    for (let i = n - 2; i >= 0; i--) {
        answer[i] = alpha[i] * answer[i + 1] + beta[i];
    }
    return answer;
}
function explicitSchema(initialTemperature, boundaryCondition) {
    const bc = boundaryCondition;
    console.log('bc: ', bc);
    const inT = [...initialTemperature];
    console.log('inT: ', inT);
    const { a, k1, k2, deltaT, deltaX } = abbr(bc);
    const r = (i) => bc.R1 + (i - 1) * deltaX;
    //init resulted vector with unphysical
    const finalT = new Array(inT.length + 2);
    //first value adding in array
    const T00 = inT[1] + 2 * k1 * deltaX * (bc.T1 - inT[0]); //(inT[0] + k1 * deltaX * bc.T1) / (1 + k1 * deltaX)
    inT.unshift(T00);
    //last value adding in array
    const T0n = inT[inT.length - 2] + 2 * k2 * deltaX * (bc.T2 - inT[inT.length - 1]); //(inT[inT.length - 1] + k2 * deltaX * bc.T2) / (1 + k2 * deltaX)
    inT.push(T0n);
    console.log('inT after adding: ', inT);
    //finding in range of [1, n-1]
    for (let i = 1; i < inT.length - 1; i++) {
        finalT[i] = inT[i] + a * deltaT / (deltaX ^ 2) * (inT[i + 1] - 2 * inT[i] + inT[i - 1]) +
            a * deltaT / (2 * deltaX * r(i)) * (inT[i + 1] - inT[i - 1]);
    }
    return finalT.slice(1, -1); //drop first and last empty array elements
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
function implicitSchema(initialTemperature, boundaryCondition) {
    const bc = boundaryCondition;
    console.log('bc: ', bc);
    const inT = [...initialTemperature];
    console.log('inT: ', inT);
    const { a, k1, k2, deltaT, deltaX } = abbr(bc);
    const r = (i) => bc.R1 + (i - 1) * deltaX;
    //init resulted vector with unphysical
    let finalT = new Array(inT.length);
    //prepare vectors of coefficients
    let A = new Array(inT.length - 1);
    let B = new Array(inT.length - 1);
    let C = new Array(inT.length - 1);
    let D = new Array(inT.length - 1);
    //first row
    A[0] = 0; //don't care
    B[0] = (-1 / (a * deltaT) - 2 / (deltaX ^ 2)) - k1 * (2 / deltaX - 1 / r(1));
    C[0] = (2 / (deltaX ^ 2));
    D[0] = inT[0] / (a * deltaT) - bc.T1 * k1 * (2 / deltaX - 1 / r(1)); //first physical initial temperature
    //i-th row
    for (let i = 1; i < (inT.length - 1); i++) {
        A[i] = 1 / (deltaX ^ 2) - 1 / (2 * deltaX * r(i));
        B[i] = -1 / (a * deltaT) - 2 / (deltaX ^ 2);
        C[i] = 1 / (deltaX ^ 2) + 1 / (2 * deltaX * r(i));
        D[i] = inT[i - 1] / (a * deltaT);
    }
    //last row
    A[inT.length - 1] = (-1 / (a * deltaT) - 2 / (deltaX ^ 2)) - k2 * (2 / deltaX + 1 / r(1));
    B[inT.length - 1] = (2 / (deltaX ^ 2));
    C[inT.length - 1] = 0; //don't care
    D[inT.length - 1] = inT[inT.length - 2] / (a * deltaT) - bc.T2 * k2 * (2 / deltaX + 1 / r(inT.length - 1)); //last physical initial temperature
    console.log("C = ", C);
    finalT = swipeMethod1(A, B, C, D);
    return finalT; //drop first and last empty array elements
}
//исполняющий код
//let temp = makeMatrixFromBC(2, boundaryCondition)
// temp[0][1][1] = 1
// console.log('temp[0][1][2]: ', temp[0][1][1]);
// console.log(temp);
console.log("result: ", explicitSchema(boundaryCondition.T0, boundaryCondition));
console.log("result: ", implicitSchema(boundaryCondition.T0, boundaryCondition));
// console.log(swipeMethod([[1, 3, 0, 0], [2, -5, 7, 0], [0, 4, 5, 8], [0, 0, 0, 6]], [1, 2, 3, 0]))
// console.log(swipeMethod1([0,2,4,0],[1,-5,5,6],[3,7,8,1],[1,2,3,0]))//correct!
