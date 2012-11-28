var Constants = {
    //Число пи
    Pi: 3.141592654,
    Pi2: 6.283185307,
    PiDiv2: 1.570796327,
    PiDiv4: 0.785398163,

    //Константы, связанные с числами с фиксированной точкой
    FIXP16_SHIFT: 16,
    FIXP16_MAG: 65536,
    FIXP16_DP_MASK: 0x0000ffff,
    FIXP16_WP_MASK: 0xffff0000,
    FIXP16_ROUND_UP: 0x00008000,

    //Малые числа (сравнения)
    Epsilon_E4: 1E-4,
    Epsilon_E5: 1E-5,
    Epsilon_E6: 1E-6,

    //Параметры пересечения отрезков
    ParamLineNoIntersect: 0,
    ParamLineIntersectInSegment: 1,
    ParamLineIntersectOutSegment: 2,
    ParamLineIntersectEverywhere: 3,

    //Единичные матрицы
    Matrix4x4Imat: {
        M: [[1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]]
    },

    Matrix4x3ImatMathIncorrect: {
        M: [[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [0, 0, 0]]
    },

    Matrix3x3Imat: {
        M: [[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]]
    },

    Matrix2x2Imat: {
        M: [[1, 0],
            [0, 1]]
    }
}