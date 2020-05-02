var v1 = [1,2,5,5,7,9];

var v2 = [1,2,4,6,7,8];


var sum = 0.0;


for (var i = 0; i < v1.length; i++){
    sum += Math.pow(v2[i] - v1[i], 2);
}

console.log(sum);
console.log(Math.pow(sum, 0.5));