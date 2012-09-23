

///////////////////////////////////////////Точки////////////////////////////////////////

//Двумерная точка
// var Point2D = new Object();
// Point2D = {x:0.0, y:0.0}
function Point2D(x, y)
{
	this.x = x;
	this.y = y;
}

function Vector2D(point2D_x, point2D_y)
{
	this.x = point2D_x;
	this.y = point2D_y;
}

// function Vector2D(point2D_x, Point2D_y)
// {
// 	this.start = point2D_x;
// 	this.end = point2D_y;
// }

//Трехмерная точка
// var Point3D = new Object();
// Point3D = {x:0.0, y:0.0, x:0.0}
function Point3D(x, y, z)
{
	this.x = x;
	this.y = y;
	this.z = z;
}

function Vector3D(x, y, z)
{
	this.x = x;
	this.y = y;
	this.z = z;
}

//Четырехмерная точка
var Point4D = new Object();
Point4D = {x:0.0, y:0.0, z:0.0, w:0.0}
function Point4D(x, y, z, w)
{
	this.x = x;
	this.y = y;
	this.z = z;
	this.w = w;
}

function Vector4D(x, y, z, w)
{
	this.x = x;
	this.y = y;
	this.z = z;
	this.w = w;
}

//////////////////////////////////////////Параметризированные прямые/////////////////////

//Двумерная параметричкская прямая
var ParamLine2D = new Object();
ParamLine2D.p0 =  Point2D; // Начальная точка
ParamLine2D.p1 =  Point2D; // Конечная точка
ParamLine2D.Vector2DasPoint2D = Point2D; // Вектор направления |v|=|p0 -> p1|
function ParamLine2D(p0, p1, vector2DasPoint2D)
{
	this.p0 = p0; // Начальная точка
	this.p1 = p1; // Конечная точка
	this.Vector2DasPoint2D = vector2DasPoint2D; // Вектор направления |v|=|p0 -> p1|
}

//Трехмерная параметрическая прямая
var ParamLine3D = new Object();
ParamLine3D.p0 = Point3D; //Начальная точка
ParamLine3D.p1 = Point3D; //Конечная точка
ParamLine3D.Vector3DasPoint3D = Point3D; //Вектор направления |v|=|p0 -> p1|
function ParamLine3D(p0, p1, vector3DasPoint3D)
{
	this.p0 = p0; // Начальная точка
	this.p1 = p1; // Конечная точка
	this.Vector3DasPoint3D = vector3DasPoint3D; //Вектор направления |v|=|p0 -> p1|
}

/////////////////////////////////////////Трехмерные плоскости/////////////////////////////

//Трехмерная плоскость
var Plane3D = new Object();
Plane3D.p0 = Point3D; // Точка на полоскости
Plane3D.NormalVectorAsPoint3D = Point3D; // Нормальный вектор
function Plane3D(p0, normalVectorAsPoint3D)
{
	this.p0 = p0; // Точка на полоскости
	this.NormalVectorAsPoint3D = normalVectorAsPoint3D; // Нормальный вектор
}

/////////////////////////////////////////Матрицы//////////////////////////////////////////

//Матрица 4х4
// var Matrix4x4 = new Object();
// Matrix4x4.M = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]];
function Matrix4x4(){
	this.M = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]];
}

//Матрица 4х3
// var Matrix4x3 = new Object();
// Matrix4x3.M = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
function Matrix4x3(){
	this.M = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
}

//Матрица 1х4
// var Matrix1x4 = new Object();
// Matrix1x4.M = [0.0, 0.0, 0.0, 0.0];
function Matrix1x4(){
	this.M = [0.0, 0.0, 0.0, 0.0];
}

//Матрица 3х3
// var Matrix3x3 = new Object();
// Matrix3x3.M = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
function Matrix3x3(){
	this.M = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
}

//Матрица 1х3
// var Matrix1x3 = new Object();
// Matrix1x3.M = [0.0, 0.0, 0.0];
function Matrix1x3(){
	this.M = [0.0, 0.0, 0.0];
}

//Матрица 3х2
//var Matrix3x2 = new Object();
//Matrix3x2.M = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];
function Matrix3x2(){
	this.M = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]];s
}

//Матрица 2х2
// var Matrix2x2 = new Object();
// Matrix2x2.M = [[0.0, 0.0], [0.0, 0.0]];
function Matrix2x2(){
	this.M = [[0.0, 0.0], [0.0, 0.0]];
}

//Матрица 1х2
// var Matrix1x2 = new Object();
// Matrix1x2.M = [0.0, 0.0];
function Matrix1x2(){
	this.M = [0.0, 0.0];
}

////////////////////////////////////////Квантернионы/////////////////////////////////////////

//Кватернион
// var QUAT = new Object();
// QUAT.M = [0.0, 0.0, 0.0, 0.0]; // Векторная часть xi + yj + zk
// QUAT.q0 = 0.0; // Действительная часть
// QUAT = {w: 0.0, x: 0.0, y:0.0, z:0.0}

function QUAT(){
	this.M = [0.0, 0.0, 0.0, 0.0]; // Векторная часть xi + yj + zk
	this.q0 = 0.0; // Действительная часть
}

///////////////////////////////////////Угловые системы координат//////////////////////////////

//Двумерные полярные координаты
var Polar2D = new Object();
Polar2D.r = 0.0; //Полярный радиус
Polar2D.theta = 0.0; //Полярный угол
function Polar2D(){
	this.r = 0.0; // Полярный радиус
	this.theta = 0.0; // Полярный угол
}

//Трехмерные цилиндрические координаты
var Cylindrical3D = new Object();
Cylindrical3D.r = 0.0; //Полярный радиус точки
Cylindrical3D.theta = 0.0; //Полярный угол точки
Cylindrical3D.z = 0.0; //z-координата точки
function Cylindrical3D(){
	this.r = 0.0; //Полярный радиус точки
	this.theta = 0.0; //Полярный угол точки
	this.z = 0.0; //z-координата точки
}

//Трехмерные сферичкские координаты
var Spherical3D = new Object();
Spherical3D.p = 0.0; //Полярный радиус
Spherical3D.theta =  0.0; //Широта
Spherical3D.phi = 0.0; //Долгота
function Spherical3D(){
	this.p = 0.0; //Полярный радиус
	this.theta =  0.0; //Широта
	this.phi = 0.0; //Долгота
}














