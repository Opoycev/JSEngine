var Constants = new Object();

//Число пи
Constants.Pi = 3.141592654;
Constants.Pi2 = 6.283185307;
Constants.PiDiv2 = 1.570796327;
Constants.PiDiv4 = 0.785398163;

//Константы, связанные с числами с фиксированной точкой
Constants.FIXP16_SHIFT = 16;
Constants.FIXP16_MAG = 65536;
Constants.FIXP16_DP_MASK = 0x0000ffff;
Constants.FIXP16_WP_MASK = 0xffff0000;
Constants.FIXP16_ROUND_UP = 0x00008000;

//Малые числа (сравнения)
Constants.Epsilon_E4 = 1E-4;
Constants.Epsilon_E5 = 1E-5;
Constants.Epsilon_E6 = 1E-6;

//Параметры пересечения отрезков
Constants.ParamLineNoIntersect = 0;
Constants.ParamLineIntersectInSegment = 1;
Constants.ParamLineIntersectOutSegment = 2;
Constants.ParamLineIntersectEverywhere = 3;

//Единичные матрицы
Constants.Matrix4x4Imat = new Object();
Constants.Matrix4x4Imat.M = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];

Constants.Matrix4x3ImatMathIncorrect = new Object();
Constants.Matrix4x3ImatMathIncorrect.M = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]];


Constants.Matrix3x3Imat = new Object();
Constants.Matrix3x3Imat.M = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

Constants.Matrix2x2Imat = new Object();
Constants.Matrix2x2Imat.M = [[1, 0], [0, 1]];





