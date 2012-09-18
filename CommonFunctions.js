//Общие функции
function Min(a, b)
{
	if (a < b) return a; else return b;
}
function Max(a, b)
{
	if (a > b) return a; else return b;
}

//Функция обмена
function Swap(a, b)
{
	 var temp;
	 temp = a; a = b; b = temp;
}

//Преобразование из радиан в градусы
function RadToDeg(rad)
{
	return rad*180/Constants.Pi;
}
//Преобразование из градусов в радианы
function DegToRad(deg)
{
	return deg*Constants.Pi/180;
}

//Случайное число из промежутка
function RandomRange(a, b)
{
	return (a + (Math.random()*1000%(b - a + 1)));
}

///////////////////////////////////////////////////////////Точки и векторы//////////////////////////////////////////////

//Обнуление координат точки/вектора
function Vector2D_Zero(vector2D)
{
	vector2D.x = 0.0;
	vector2D.y = 0.0;
	return (vector2D);
}
function Vector3D_Zero(vector3D)
{
	vector3D.x = 0.0;
	vector3D.y = 0.0;
	vector3D.z = 0.0;
	return (vector3D);
}
function Vector4D_Zero(vector4D)
{
	vector4D.x = 0.0;
	vector4D.y = 0.0;
	vector4D.z = 0.0;
	vector4D.w = 1.0;
}

//Инициализация векторов
function Vector2D_InitXY(vector2D, x, y)
{
	vector2D.x = x;
	vector2D.y = y;
	return vector2D;
}
function Vector3D_InitXYZ(vector3D, x, y, z)
{
	vector3D.x = x;
	vector3D.y = y;
	vector3D.z = z;
	return vector3D;
}
function Vector4D_InitXYZ(vector4D, x, y, z)
{
	vector4D.x = x;
	vector4D.y = y;
	vector4D.z = z;
	vector4D.w = 1.0;
	return vector4D;
}

//Инициализация векторов другими векторами
function Vector2D_Init(vector2Ddest, vector2Dsource)
{
	vector2Ddest.x = vector2Dsource.x;
	vector2Ddest.y = vector2Dsource.y;
	return vector2Ddest;
}
function Vector3D_Init(vector3Ddest, vector3Dsource)
{
	vector3Ddest.x = vector3Dsource.x;
	vector3Ddest.y = vector3Dsource.y;
	vector3Ddest.z = vector3Dsource.z;
	return vector3Ddest;
}
function Vector4D_Init(vector4Ddest, vector4Dsource)
{
	vector4Ddest.x = vector4Dsource.x;
	vector4Ddest.y = vector4Dsource.y;
	vector4Ddest.z = vector4Dsource.z;
	vector4Ddest.w = vector4Dsource.w;
	return vector4Ddest;
}

//Копирование векторов
function Vector2D_Copy(vector2Ddest, vector2Dsource)
{
	vector2Ddest.x = vector2Dsource.x;
	vector2Ddest.y = vector2Dsource.y;
}
function Vector3D_Copy(vector3Ddest, vector3Dsource)
{
	vector3Ddest.x = vector3Dsource.x;
	vector3Ddest.y = vector3Dsource.y;
	vector3Ddest.z = vector3Dsource.z;
}
function Vector4D_Copy(vector4Ddest, vector4Dsource)
{
	vector4Ddest.x = vector4Dsource.x;
	vector4Ddest.y = vector4Dsource.x;
	vector4Ddest.z = vector4Dsource.z;
	vector4Ddest.w = vector4Dsource.w;
}

//Инициализация точек
function Point2D_Init(point2Ddest, point2Dsource)
{
	point2Ddest.x = point2Dsource.x;
	point2Ddest.y = point2Dsource.y;
}
function Point3D_Init(point3Ddest, point3Dsource)
{
	point3Ddest.x = point3Dsource.x;
	point3Ddest.y = point3Dsource.y;
	point3Ddest.z = point3Dsource.z;
}
function Point4D_Init(point4Ddest, point4Dsource)
{
	point4Ddest.x = point4Dsource.x;
	point4Ddest.y = point4Dsource.y;
	point4Ddest.z = point4Dsource.z;
	point4Ddest.w = point4Dsource.w;
}

//Копирование точек
function Point2D_Copy(point2Ddest, point2Dsource)
{
	point2Ddest.x = point2Dsource.x;
	point2Ddest.y = point2Dsource.y;
}
function Point3D_Copy(point3Ddest, point3Dsource)
{
	point3Ddest.x = point3Dsource.x;
	point3Ddest.y = point3Dsource.y;
	point3Ddest.z = point3Dsource.z;
}
function Point4D_Copy(point4Ddest, point4Dsource)
{
	point4Ddest.x = point4Dsource.x;
	point4Ddest.y = point4Dsource.y;
	point4Ddest.z = point4Dsource.z;
	point4Ddest.w = point4Dsource.w;
}

///////////////////////////////////////////////////////Матричные функции////////////////////////////////////////////////

//Функции обнуления матриц
function Matrix_Zero_2x2(matrix2x2)
{
	matrix2x2.M[0][0] = 0.0; matrix2x2.M[0][1] = 0.0;
	matrix2x2.M[1][0] = 0.0; matrix2x2.M[1][1] = 0.0;
}
function Matrix_Zero_3x3(matrix3x3)
{
	matrix3x3.M[0][0] = 0.0; matrix3x3.M[0][1] = 0.0; matrix3x3.M[0][2] = 0.0;
	matrix3x3.M[1][0] = 0.0; matrix3x3.M[1][1] = 0.0; matrix3x3.M[1][2] = 0.0;
	matrix3x3.M[2][0] = 0.0; matrix3x3.M[2][1] = 0.0; matrix3x3.M[2][2] = 0.0;
}
function Matrix_Zero_4x4(matrix4x4)
{
	matrix4x4.M[0][0] = 0.0; matrix4x4.M[0][1] = 0.0; matrix4x4.M[0][2] = 0.0; matrix4x4.M[0][3] = 0.0;
	matrix4x4.M[1][0] = 0.0; matrix4x4.M[1][1] = 0.0; matrix4x4.M[1][2] = 0.0; matrix4x4.M[1][3] = 0.0;
	matrix4x4.M[2][0] = 0.0; matrix4x4.M[2][1] = 0.0; matrix4x4.M[2][2] = 0.0; matrix4x4.M[2][3] = 0.0;
	matrix4x4.M[3][0] = 0.0; matrix4x4.M[3][1] = 0.0; matrix4x4.M[3][2] = 0.0; matrix4x4.M[3][3] = 0.0;
}
function Matrix_Zero_4x3(matrix4x3)
{
	matrix2x2.M[0][0] = 0.0; matrix2x2.M[0][1] = 0.0; matrix2x2.M[0][2] = 0.0;
	matrix2x2.M[1][0] = 0.0; matrix2x2.M[1][1] = 0.0; matrix2x2.M[1][2] = 0.0;
	matrix2x2.M[2][0] = 0.0; matrix2x2.M[2][1] = 0.0; matrix2x2.M[2][2] = 0.0;
	matrix2x2.M[3][0] = 0.0; matrix2x2.M[3][1] = 0.0; matrix2x2.M[3][2] = 0.0;
}

//Функции инициализации единичными матрицами
function Matrix_Identity_2x2(matrix2x2)
{
	matrix2x2 = Constants.Matrix2x2Imat;
}
function Matrix_Identity_3x3(matrix3x3)
{
	matrix3x3 = Constants.Matrix3x3Imat;
}
function Matrix_Identity_4x4(matrix4x4)
{
	matrix4x4 = Constants.Matrix4x4Imat;
}
function Matrix_Identity_4x3(matrix4x3)
{
	matrix4x3 = Constants.Matrix4x3ImatMathIncorrect;
}

//Функции копирования матриц
function Matrix_Copy_2x2(matrix2x2dest, matrix2x2source)
{
	matrix2x2dest.M[0][0] = matrix2x2source.M[0][0]; matrix2x2dest.M[0][1] = matrix2x2source.M[0][1];
	matrix2x2dest.M[1][0] = matrix2x2source.M[1][0]; matrix2x2dest.M[1][1] = matrix2x2source.M[1][1];
}
function Matrix_Copy_3x3(matrix3x3dest, matrix3x3source)
{
	matrix2x2dest.M[0][0] = matrix2x2source.M[0][0]; matrix2x2dest.M[0][1] = matrix2x2source.M[0][1]; matrix2x2dest.M[0][2] = matrix2x2source.M[0][2];
	matrix2x2dest.M[1][0] = matrix2x2source.M[1][0]; matrix2x2dest.M[1][1] = matrix2x2source.M[1][1]; matrix2x2dest.M[1][2] = matrix2x2source.M[1][2];
	matrix2x2dest.M[2][0] = matrix2x2source.M[2][0]; matrix2x2dest.M[2][1] = matrix2x2source.M[2][1]; matrix2x2dest.M[2][2] = matrix2x2source.M[2][2];
}
function Matrix_Copy_4x4(matrix4x4dest, matrix4x4source)
{
	matrix2x2dest.M[0][0] = matrix2x2source.M[0][0]; matrix2x2dest.M[0][1] = matrix2x2source.M[0][1]; matrix2x2dest.M[0][2] = matrix2x2source.M[0][2]; matrix2x2dest.M[0][3] = matrix2x2source.M[0][3]; 
	matrix2x2dest.M[1][0] = matrix2x2source.M[1][0]; matrix2x2dest.M[1][1] = matrix2x2source.M[1][1]; matrix2x2dest.M[1][2] = matrix2x2source.M[1][2]; matrix2x2dest.M[1][3] = matrix2x2source.M[1][3];
	matrix2x2dest.M[2][0] = matrix2x2source.M[2][0]; matrix2x2dest.M[2][1] = matrix2x2source.M[2][1]; matrix2x2dest.M[2][2] = matrix2x2source.M[2][2]; matrix2x2dest.M[2][3] = matrix2x2source.M[2][3]; 
	matrix2x2dest.M[3][0] = matrix2x2source.M[3][0]; matrix2x2dest.M[3][1] = matrix2x2source.M[3][1]; matrix2x2dest.M[3][2] = matrix2x2source.M[3][2]; matrix2x2dest.M[3][3] = matrix2x2source.M[3][3];  
}
function Matrix_Copy_4x3(matrix4x3dest, matrix4x3source)
{
	matrix2x2dest.M[0][0] = matrix2x2source.M[0][0]; matrix2x2dest.M[0][1] = matrix2x2source.M[0][1]; matrix2x2dest.M[0][2] = matrix2x2source.M[0][2];
	matrix2x2dest.M[1][0] = matrix2x2source.M[1][0]; matrix2x2dest.M[1][1] = matrix2x2source.M[1][1]; matrix2x2dest.M[1][2] = matrix2x2source.M[1][2];
	matrix2x2dest.M[2][0] = matrix2x2source.M[2][0]; matrix2x2dest.M[2][1] = matrix2x2source.M[2][1]; matrix2x2dest.M[2][2] = matrix2x2source.M[2][2];
	matrix2x2dest.M[3][0] = matrix2x2source.M[3][0]; matrix2x2dest.M[3][1] = matrix2x2source.M[3][1]; matrix2x2dest.M[3][2] = matrix2x2source.M[3][2];
}

//Транспонирование матриц
function MAT_TRANSPOSE_3x3(matrix3x3)
{
	var mt = Matrix3x3;
	mt.M[0][0] = matrix3x3.M[0][0]; mt.M[0][1] = matrix3x3.M[1][0]; mt.M[0][2] = matrix3x3.M[2][0];
	mt.M[1][0] = matrix3x3.M[0][1]; mt.M[1][1] = matrix3x3.M[1][1]; mt.M[1][2] = matrix3x3.M[2][1];
	mt.M[2][0] = matrix3x3.M[0][2]; mt.M[2][1] = matrix3x3.M[1][2]; mt.M[2][2] = matrix3x3.M[2][2];
	
	return mt;
}
function MAT_TRANSPOSE_4x4(matrix4x4)
{
	var mt = Matrix4x4;
	mt.M[0][0] = matrix4x4.M[0][0]; mt.M[0][1] = matrix4x4.M[1][0];
	mt.M[0][2] = matrix4x4.M[2][0]; mt.M[0][3] = matrix4x4.M[3][0]; 
	mt.M[1][0] = matrix4x4.M[0][1]; mt.M[1][1] = matrix4x4.M[1][1]; 
	mt.M[1][2] = matrix4x4.M[2][1]; mt.M[1][3] = matrix4x4.M[3][1]; 
	mt.M[2][0] = matrix4x4.M[0][2]; mt.M[2][1] = matrix4x4.M[1][2]; 
	mt.M[2][2] = matrix4x4.M[2][2]; mt.M[2][3] = matrix4x4.M[3][2]; 
	mt.M[3][0] = matrix4x4.M[0][3]; mt.M[3][1] = matrix4x4.M[1][3]; 
	mt.M[3][2] = matrix4x4.M[2][3]; mt.M[3][3] = matrix4x4.M[3][3];
	matrix4x4 = mt;   
}
function MAT_TRANSPOSE_3X3_FROM_MATRIX(matrix3x3dest, matrix3x3source)
{
	matrix3x3dest.M[0][0] = matrix3x3source.M[0][0]; matrix3x3dest.M[0][1] = matrix3x3source.M[1][0]; matrix3x3dest.M[0][2] = matrix3x3source.M[2][0];
	matrix3x3dest.M[1][0] = matrix3x3source.M[0][1]; matrix3x3dest.M[1][1] = matrix3x3source.M[1][1]; matrix3x3dest.M[1][2] = matrix3x3source.M[2][1];
	matrix3x3dest.M[2][0] = matrix3x3source.M[0][2]; matrix3x3dest.M[2][1] = matrix3x3source.M[1][2]; matrix3x3dest.M[2][2] = matrix3x3source.M[2][2];
}
function MAT_TRANSPOSE_4X4_FROM_MATRIX(matrix4x4dest, matrix4x4source)
{
	matrix4x4dest.M[0][0] = matrix4x4source.M[0][0]; matrix4x4dest.M[0][1] = matrix4x4source.M[1][0]; 
	matrix4x4dest.M[0][2] = matrix4x4source.M[2][0]; matrix4x4dest.M[0][3] = matrix4x4source.M[3][0]; 
	matrix4x4dest.M[1][0] = matrix4x4source.M[0][1]; matrix4x4dest.M[1][1] = matrix4x4source.M[1][1];
	matrix4x4dest.M[1][2] = matrix4x4source.M[2][1]; matrix4x4dest.M[1][3] = matrix4x4source.M[3][1]; 
	matrix4x4dest.M[2][0] = matrix4x4source.M[0][2]; matrix4x4dest.M[2][1] = matrix4x4source.M[1][2]; 
	matrix4x4dest.M[2][2] = matrix4x4source.M[2][2]; matrix4x4dest.M[2][3] = matrix4x4source.M[3][2]; 
	matrix4x4dest.M[3][0] = matrix4x4source.M[0][3]; matrix4x4dest.M[3][1] = matrix4x4source.M[1][3]; 
	matrix4x4dest.M[3][2] = matrix4x4source.M[2][3]; matrix4x4dest.M[3][3] = matrix4x4source.M[3][3]; 
}

//Обмен столбцов
// [Tested]
function MAT_COLUMN_SWAP_2X2(matrix2x2, c, matrix1x2)
{
	matrix2x2.M[0][c] = matrix1x2.M[0]; matrix2x2.M[1][c] = matrix1x2.M[1];
}
function MAT_COLUMN_SWAP_3X3(matrix3x3, c, matrix1x3)
{
	matrix3x3.M[0][c] = matrix1x3.M[0]; matrix3x3.M[1][c] = matrix1x3.M[1];
	matrix3x3.M[2][c] = matrix1x3.M[2];
}
function MAT_COLUMN_SWAP_4X4(matrix4x4, c, matrix1x4)
{
	matrix4x4.M[0][c] = matrix1x4.M[0]; matrix4x4.M[1][c] = matrix1x4.M[1];
	matrix4x4.M[2][c] = matrix1x4.M[2]; matrix4x4.M[3][c] = matrix1x4.M[3];
}
function MAT_COLUMN_SWAP_4X3(matrix4x3, c, matrix1x4)
{
	matrix4x3.M[0][c] = matrix1x4.M[0]; matrix4x3.M[1][c] = matrix1x4.M[1];
	matrix4x3.M[2][c] = matrix1x4.M[2]; matrix4x3.M[3][c] = matrix1x4.M[3];
}

/////////////////////////////////////////////////////////////Кватернионы/////////////////////////////////////////////////////////

function QUAT_ZERO(quat)
{
	quat.w = 0.0; quat.M[0] = 0.0;
	quat.x = 0.0; quat.M[1] = 0.0;
	quat.y = 0.0; quat.M[2] = 0.0;
	quat.z = 0.0; quat.M[3] = 0.0;
}
function QUAT_InitWXYZ(quat, w, x, y, z)
{
	quat.w = w; quat.M[0] = w;
	quat.x = x; quat.M[1] = x;
	quat.y = y; quat.M[2] = y;
	quat.z = z; quat.M[3] = z;
}
function QUAT_Init_Vector3D(quat, vector)
{
	quat.w = 0.0; quat.x = vector.x;
	quat.y = vector.y; quat.z = vector.z;
}
function QUAT_Init(quatDest, quatSource)
{
	quatDest.w = quatSource.w; quatDest.x = quatSource.x;
	quatDest.y = quatSource.y; quatDest.z = quatSource.z;
}
function QUAT_Copy(quatDest, quatSource)
{
	quatDest.w = quatSource.w; quatDest.x = quatSource.x;
	quatDest.y = quatSource.y; quatDest.z = quatSource.z;
}

///////////////////////////////////////////////////////Математика с фиксированной точкой////////////////////////////////////////////
//!!!!!!!Под вопросом


//Выделение целой и дробной частей числа с фиксированной точкой в формате 16.16
function FIX16_WP(fp)
{
	return (fp >> Constants.FIXP16_SHIFT);
}
function FIX16_DP(fp)
{
	return (fp && Constants.FIXP16_DP_MASK);
}
// Преобразование целых чисел и чисел с плавающей точкой в числа с фиксированной точкой в формате 16.16
function INT_TO_FIX16(i)
{
	return (i << Constants.FIXP16_SHIFT);
}
function FLOAT_TO_FIX16(f)
{
	return (f * Constants.FIXP16_MAG + 0.5);
}
// Преобразование числа с фиксированной точкой в число с плавающей точкой
function FIXP16_TO_FLOAT(f)
{
	return (f / Constants.FIXP16_MAG);
}




























































