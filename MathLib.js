MathLib = {
	FastDistance2D: function(x, y)
	{
		//вычисляет расстояние между точками 0,0 и x, y с погрешностью 3,5%
		
		//Вычисляем абсолютные значения
		var x = Math.abs(x);
		var y = Math.abs(y);
		
		var mn = Min(x, y);
		
		return (x + y - (mn >> 1) - (mn >> 2) + (mn >> 4));
			
	},
	
	FastDistance3D: function(fx, fy, fz)
	{
		var x, y, z;
		var temp;
		
		x = Math.abs(fx) * 1024;
		y = Math.abs(fy) * 1024;
		z = Math.abs(fz) * 1024;
		
		if (y < x) {temp = y; y = x; x = temp;}; //Swap(x, y, temp);
		if (z < y) {temp = z; z = y; y = temp;}; //Swap(y, z, temp);
		if (y < x) {temp = y; y = x; x = temp;}; //Swap(x, y, temp);
		
		var dist = (z + 11 * (y >> 5) + (x >> 2));
		
		return (dist >> 10);
	},
	
	Polar2DToPoint2D: function(polar, rect)
	{
		//Преобразование из полярных в декартову систему координат
		rect.x = polar.r * Math.cos(polar.theta);
		rect.y = polar.r * Math.sin(polar.theta);
	},
	
	Polar2D_to_RectXY: function(polar, x, y)
	{
		//Преобразование из полярных координат в декартовы
		x = polar.r * Math.cos(polar.theta);
		y = polar.r * Math.sin(polar.theta);
	},
	
	Point2D_to_Polar2D: function(rect, polar)
	{
		//Преобразование из двумерных декартовых координат в полярные
		polar.r = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		polar.theta = Math.atan(rect.y / rect.x);
	},
	
	POINT2D_To_PolarRTr: function(rect, r, theta)
	{
		//Функция преобразует двумерные декартовы координаты в явные значения полярных координат.
		r = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		theta = Math.atan(rect.y / rect.x);
	},
	
	CYLINDRICAL3D_To_POINT3D: function(cyl, rect)
	{
		//Данная функция преобразует цилиндрические координаты в декартовы.
		rect.x = cyl.r * Math.cos(cyl.theta);
		rect.y = cyl.r * Math.sin(cyl.theta);
		rect.z = cyl.z;
	},
	
	CYLINDRICAL3D_To_RectXYZ: function(cyl, x, y, z)
	{
		//Данная функция преобразует цилиндрические координаты в явные значения декартовых координат.
		x = cyl.r * Math.cos(cyl.theta);
		y = cyl.r * Math.sin(cyl.theta);
		z = cyl.z;
	},
	
	POINT3D_To_CYLINDRICAL3D: function(rect, cyl)
	{
		//Функция преобразует декартовы координаты в цилиндрические.
		cyl.r = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		cyl.theta = Math.atan(rect.y / rect.x);
		cyl.z = rect.z;
	},
	
	POINT3D_To_CylindricalRThZ: function(rect, r, theta, z)
	{
		//Функция преобразует декартовы координаты в явные значения цилиндрических координат
		r = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		theta = Math.atan(rect.y / rect.x);
		z = rect.z;
	},
	
	SPHERICAL3D_ToPDINT3D: function(sph, rect)
	{
		//Функция преобразует сферические координаты в декартовы
		var r = 0.0;
		r = sph.p * Math.sin(sph.phi);
		rect.z = sph.p * Math.cos(sph.phi);
		
		rect.x = r * Math.cos(sph.theta);
		rect.y = r * Math.sin(sph.theta);
	},
	
	SPHERICAL3D_To_RectXYZ: function(sph, x, y, z)
	{
		// Функция преобразует сферические координаты в явные значения декартовых координат
		var r = 0.0;
		r = sph.p * Math.sin(sph.phi);
		z = sph.p * Math.sin(sph.thi);
		
		x = r * Math.cos(sph.theta);
		y = r * Math.sin(sph.theta);
	},
	
	POINT3D_To_SPHERICAL3D: function(rect, sph)
	{
		sph.p = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		sph.theta = Math.atan(rect.y / rect.x);
		
		var r = sph.p * Math.sin(sph.phi);
		sph.phi = Math.asin(r / sph.p);
	},
	
	POINT3D_To_SphericalPThPh: function(rect, p, theta, phi)
	{
		//Функция преобразует декартовы координаты в явные значения сферических координат
		p = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y) + (rect.z * rect.z));
		theta = Math.atan(rect.y / rect.x);
		
		var r = Math.sqrt((rect.x * rect.x) + (rect.y * rect.y));
		phi = Math.asin(r / p);
	},
	
	
	////////////////////////////////////////Функции для работы с векторами//////////////////////////////////////
	//Функции VECTOR*D_Add() суммируют переданные в качестве параметров векторы и возвращают полученный в результате вектор
	VECTOR2D_Add: function(va, vb, vsum)
	{
		vsum.x = va.x + vb.x;
		vsum.y = va.y + vb.y; 
	},
	
	VECTOR2D_AddAndReturnVector2D: function(va, vb)
	{
		var vsum = new Vector2D();

		vsum.x = va.x + vb.x;
		vsum.y = va.y + vb.y;
	
		return vsum;
	},
	
	VECTOR2D_Sub: function(va, vb, vdiff)
	{
		vdiff.x = va.x - vb.x;
		vdiff.y = vb.y - vb.y;
	},
	
	VECTOR2D_SubAndReturnVector2D: function(va, vb)
	{
		var vdiff = new Vector2D(0, 0);
		
		vdiff.x = va.x - vb.x;
		vdiff.y = va.y - vb.y;
		
		return vdiff;
	},
	
	VECTOR2D_Scale: function(k, va, vscaled)
	{
		//Функции void VECTOR*D_Scale() масштабируют переданный в качестве параметра va вектор, увеличивая его в k раз, и возвращают масштабированный вектор как vscaled.
		vscaled.x = k * va.x;
		vscaled.y = k * va.y;
	},
	
	VECTOR2D_Scale_Single: function(k, va)
	{
		va.x *= k;
		va.y *= k; 
	},
	
	VECTOR2D_Dot: function(va, vb)
	{
		//Функции float VECTOR*D_Dot() возвращают результат скалярного произведения двухвекторов. В четырехмерной версии компонента w в вычислении не участвует
		return ((va.x * vb.x) + (va.y * vb.y));
	},
	
	VECTOR2D_Length: function(va)
	{
		return(Math.sqrt(va.x * va.x) + (va.y * va.y));
	},
	
	VECTOR2D_Length_Fast: function(va)
	{
		//!!!!!!!Округление
		return(FastDistance2D(va.x, va.y));
	},
	
	 VECTOR2D_Normalize: function(va)
	 {
	 	//Нормализация вектора
	 	var length = Math.sqrt(va.x * va.x + va.y * va.y);
	 	
	 	if (length < Constants.Epsilon_E5) return;
	 	
	 	var length_env = (1.0 / length);
	 	
	 	va.x = va.x * length_inv;
	 	va.y = va.y * length_inv;
	 },
	 
	 VECTOR2D_NormalizeAndReturnInVN: function(va, vn)
	 {
	 	//Нормализация вектора и его возврат через vn
	 	VeVector2D_Zero(vn);
	 	
	 	var length = Math.sqrt(va.x * va.x + va.y * va.y);
	 	
	 	if (length < Constants.Epsilon_E5) return;
	 	
	 	var length_inv = 1.0 /length;
	 	
	 	vn.x = va.x * length_inv;
	 	vn.y = va.y * length_inv;
	 },
	 
	 VECTOR2D_Build: function(init, term, result)
	 {
	 	//Функции  VECTQR*D_Birild() строят вектор init~>term и сохраняют его в переменной result. Это хороший пример функций, которые могут работать как с векторами, так и с точками —для создания вектора, определяемого двумя точками
	 	result.x = term.x - init.x;
	 	result.y = term.y - init.y;
	 },
	 
	 VECTOR2D_CosTh: function(va, vb)
	 {
	 	//Функции  VECTOR*D_CosTh() вычисляют косинус угла между двумя векторами
	 	return(MathLib.VECTOR2D_Dot(va, vb)/(MathLib.VECTOR2D_Length(va) * MathLib.VECTOR2D_Length(vb)));
	 },
	 
	 Vector2D_Print: function(va, name, outputTag)
	 {
	 	var OutputNode = document.getElementById(outputTag);
		var val = document.createElement("p");
		var OutputString = "";
		OutputString = OutputString.concat(name, ": x=", va.x.toString(), ", y=", va.y.toString());		
		val.innerHTML = OutputString;
		OutputNode.appendChild(val);
	 },
	 
	 VECTOR3D_Add: function(va, vb, vsum)
	 {
	 	//Сложение векторов
	 	vsum.x = va.x + vb.x;
	 	vsum.y = va.y + vb.y;
	 	vsum.z = va.z + vb.z;
	 },
	 
	 VECTOR3D_AddAndReturnVector3D: function(va, vb)
	 {
	 	var vsum = new Vector3D(0, 0, 0);
	 	
	 	vsum.x = va.x + vb.x;
	 	vsum.y = va.y + vb.y;
	 	vsum.z = va.z + vb.z;
	 	
	 	return vsum;
	 },
	 
	 VECTOR3D_Sub: function(va, vb, vdiff)
	 {
	 	//Вычитание векторов
	 	vdiff.x = va.x - vb.x;
	 	vdiff.y = va.y - vb.y;
	 	vdiff.z = va.z - vb.z;
	 },
	 
	 VECTOR3D_SubAndReturnVector3D: function(va, vb)
	 {
	 	var vdiff = new Vector3D(0, 0, 0);
	 	
	 	vdiff.x = va.x - vb.x;
	 	vdiff.y = va.y - vb.y;
	 	vdiff.z = va.z - vb.z;
	 	
	 	return vdiff;
	 },
	 
	 VECTOR3D_Scale: function(k, va)
	 {
	 	//Масштабирование вектора по коэффициенту к
	 	va.x *= k;
	 	va.y *= k;
	 	va.z *= k;
	 },
	 
	 VECTOR3D_ScaleAndReturnVScaled: function(k, va, vscaled)
	 {
	 	vscaled.x = k * va.x;
	 	vscaled.y = k * va.y;
	 	vscaled.z = k * va.z;
	 },
	 
	 VECTOR3D_Dot: function(va, vb)
	 {
	 	//Скалярное произведение
	 	return((va.x * vb.x) + (va.y * vb.y) + (va.z * vb.z));
	 },
	 
	 VECTOR3D_Cross: function(va, vb, vn)
	 {
	 	//Функции void VECTOR*D_Cross() вычисляют векторное произведение переданных в качестве параметров векторов и сохраняют результат по адресу vn. В четырехмерной версии компонента w в вычислении не участвует
	 	vn.x = ((va.y * vb.z) - (va.z * vb.y));
	 	vn.y = -((va.x * vb.z) - (va.z * vb.x));
	 	vn.z = ((va.y * vb.y) - (va.y * vb.x));
	 },
	 
	 VECTOR3D_CrossAndReturnVector3D: function(va, vb)
	 {
	 	var vn = new Vector3D(0, 0, 0);
	 	
	 	vn.x =  ( (va.y * vb.z) - (va.z * vb.y) );
		vn.y = -( (va.x * vb.z) - (va.z * vb.x) );
		vn.z =  ( (va.x * vb.y) - (va.y * vb.x) ); 

		return(vn);
	 },
	 
	 VECTOR3D_Length: function(va)
	 {
	 	return(Math.sqrt(va.x * va.x) + (va.y * va.y) + (va.z * va.z));
	 },
	 
	 VECTOR3D_Normalize: function(va)
	 {
	 	var length = Math.sqrt(va.x * va.x) + (va.y * va.y) + (va.z * va.z);
	 	
	 	if (length < Constants.Epsilon_E5) return;
	 	
	 	var length_inv = 1.0 / length;
	 	
	 	va.x *= length_inv;
		va.y *= length_inv;
		va.z *= length_inv;
	 },
	 
	 VECTOR3D_NormalizeAndReturnVN: function(va, vn)
	 {
	 	Vector3D_Zero(vn);
	 	
	 	var length = VECTOR3D_Length(va);

		if (length < EPSILON_E5) return;

		var length_inv = 1.0/length;

		vn.x = va.x * length_inv;
		vn.y = va.y * length_inv;
		vn.z = va.z * length_inv;
	 },
	 
	 VECTOR3D_Build: function(init, term, result)
	 {
	 	result.x = term.x - init.x;
	 	result.y = term.y - init.y;
	 	result.z = term.z - init.x;
	 },
	 
	 VECTOR3D_CosTh: function(va, vb)
	 {
	 	// возвращает косинус угла мужду двумя векторами
	 	return( MathLib.VECTOR3D_Dot(va, vb) / (MathLib.VECTOR3D_Length(va) * MathLib.VECTOR3D_Length(vb)) );
	 },
	 
	 ///////////////////////////////////////////////////4D версии функций//////////////////////////////////////////////////////////////////////////
	 VECTOR4D_Build: function(VECTOR4D_init, VECTOR4D_term, VECTOR4D_result)
	 {
		VECTOR4D_result.x = VECTOR4D_term.x - VECTOR4D_init.x;
		VECTOR4D_result.y = VECTOR4D_term.y - VECTOR4D_init.y;
		VECTOR4D_result.z = VECTOR4D_term.z - VECTOR4D_init.z;
		VECTOR4D_result.w = 1;

	 },
	 
	 VECTOR4D_Add: function(VECTOR4D_va, VECTOR4D_vb, VECTOR4D_vsum)
	 {
		// прибавляет va к vb и возвращает через vsum
		VECTOR4D_vsum.x = VECTOR4D_va.x + VECTOR4D_vb.x;
		VECTOR4D_vsum.y = VECTOR4D_va.y + VECTOR4D_vb.y;
		VECTOR4D_vsum.z = VECTOR4D_va.z + VECTOR4D_vb.z;
		VECTOR4D_vsum.w = 1;
	 },
	 
	 VECTOR4D_AddAndReturnVector4D: function(VECTOR4D_va, VECTOR4D_vb)
	 {
		// прибавляет ya к yb и взвращает результат через стек
		var vsum = new VECTOR4D();

		vsum.x = VECTOR4D_va.x + VECTOR4D_vb.x;
		vsum.y = VECTOR4D_va.y + VECTOR4D_vb.y;
		vsum.z = VECTOR4D_va.z + VECTOR4D_vb.z;
		vsum.w = 1;

		return(vsum);

	 },
	 
	 VECTOR4D_Sub: function(VECTOR4D_va, VECTOR4D_vb, VECTOR4D_vdiff)
	 {
		// вычитает va из vb и возвращает результат через vdiff
		VECTOR4D_vdiff.x = VECTOR4D_va.x - VECTOR4D_vb.x;
		VECTOR4D_vdiff.y = VECTOR4D_va.y - VECTOR4D_vb.y;
		VECTOR4D_vdiff.z = VECTOR4D_va.z - VECTOR4D_vb.z;
		VECTOR4D_vdiff.w = 1;
	 },
	 
	 VECTOR4D_SubAndReturnVector4D: function(VECTOR4D_va, VECTOR4D_vb)
	 {
		// вычитает va из vb и возвращает результат через стек
		var vdiff = new VECTOR4D();

		vdiff.x = VECTOR4D_va.x - VECTOR4D_vb.x;
		vdiff.y = VECTOR4D_va.y - VECTOR4D_vb.y;
		vdiff.z = VECTOR4D_va.z - VECTOR4D_vb.z;
		vdiff.w = 1;

		return(vdiff);                      
	 },
	 
	 VECTOR4D_Scale: function(/*коэффициента масштабирования*/ k, VECTOR4D_va)
	 {
		// масштабирует вектор на заданный коэффициент k
		
		// просто умножаем каждый компонент к коэффициент
		VECTOR4D_va.x*=k;
		VECTOR4D_va.y*=k;
		VECTOR4D_va.z*=k;
		VECTOR4D_va.w = 1; // остается неизменным
	 },
	 
	 VECTOR4D_Scale: function(/*коэффициента масштабирования*/k, VECTOR4D_va, VECTOR4D_vscaled)
 	 {
		// масштабирует вектор на заданный коэффицент k сохраняя оригинал

		// перемножаем каждый компонент на коэффициент
		VECTOR4D_vscaled.x = k*VECTOR4D_va.x;
		VECTOR4D_vscaled.y = k*VECTOR4D_va.y;
		VECTOR4D_vscaled.z = k*VECTOR4D_va.z;
		VECTOR4D_vscaled.w = 1

	 },

	 VECTOR4D_Dot: function(VECTOR4D_va, VECTOR4D_vb)
	 {
		// вычисляет скалярное произведение va и vb
		return( (VECTOR4D_va.x * VECTOR4D_vb.x) + (VECTOR4D_va.y * VECTOR4D_vb.y) + (VECTOR4D_va.z * VECTOR4D_vb.z) );
	 },
	 
	 VECTOR4D_Cross: function(VECTOR4D_va, VECTOR4D_vb, VECTOR4D_vn)
 	{
		// вычисляет векторное произведение между va и vb
		// и возвращает результат через vn

		VECTOR4D_vn.x =  ( (VECTOR4D_va.y * VECTOR4D_vb.z) - (VECTOR4D_va.z * VECTOR4D_vb.y) );
		VECTOR4D_vn.y = -( (VECTOR4D_va.x * VECTOR4D_vb.z) - (VECTOR4D_va.z * VECTOR4D_vb.x) );
		VECTOR4D_vn.z =  ( (VECTOR4D_va.x * VECTOR4D_vb.y) - (VECTOR4D_va.y * VECTOR4D_vb.x) ); 
		VECTOR4D_vn.w = 1; // остается неизменным
	 },
	 
	 VECTOR4D_Length: function(VECTOR4D_va)
 	 {
		// вычисляет величину вектора

		return(Math.sqrt(VECTOR4D_va.x*VECTOR4D_va.x + VECTOR4D_va.y*VECTOR4D_va.y + VECTOR4D_va.z*VECTOR4D_va.z) );
	 },
	 
	 VECTOR4D_Length_Fast: function(VECTOR4D_va)
	 {
		// вычисляет величину вектора
		return( MathLib.FastDistance3D(VECTOR4D_va.x, VECTOR4D_va.y, VECTOR4D_va.z) );
	 },
	 
	 VECTOR4D_Normalize: function(VECTOR4D_va)
	 {
		// нормализует вектор

		var length = Math.sqrt(VECTOR4D_va.x*VECTOR4D_va.x + VECTOR4D_va.y*VECTOR4D_va.y + VECTOR4D_va.z*VECTOR4D_va.z);

		// есди вектор нулевой длины, то выходим
		if (length < Constants.Epsilon_E5) 
   			return;

		var length_inv = 1.0/length;

		VECTOR4D_va.x*=length_inv;
		VECTOR4D_va.y*=length_inv;
		VECTOR4D_va.z*=length_inv;
		VECTOR4D_va.w = 1;
	 },
	 
	 VECTOR4D_NormalizeAndReturnInVN: function(VECTOR4D_va, VECTOR4D_vn)
	 {
		// нормализует вектор и возвращает рузельтат через vn

		VECTOR4D_ZERO(vn);

		var length = sqrt(VECTOR4D_va.x*VECTOR4D_va.x + VECTOR4D_va.y*VECTOR4D_va.y + VECTOR4D_va.z*VECTOR4D_va.z);

		if (length < Constants.Epsilon_E5) 
   			return;

		var length_inv = 1.0/length;

		VECTOR4D_vn.x = VECTOR4D_va.x*length_inv;
		VECTOR4D_vn.y = VECTOR4D_va.y*length_inv;
		VECTOR4D_vn.z = VECTOR4D_va.z*length_inv;
		VECTOR4D_vn.w = 1;
	 },
	 
	 VECTOR4D_CosTh: function(VECTOR4D_va, VECTOR4D_vb)
	 {
		// Вычисляет косинус угла между двумя векторами
		return(MathLib.VECTOR4D_Dot(va,vb)/(MathLib.VECTOR4D_Length(va)*MathLib.VECTOR4D_Length(vb)));
	 },
	 
	 //////////////////////////////////////////Матричные функции/////////////////////////////////////////////////////
	 //////////////////////////////////////Функции для работы с двумерными матрицами/////////////////////////////////
	 
	 Mat_Init_2X2: function(MATRIX2X2_ma, m00, m01, m10, m11)
	 {
		// заполняет матрицу 2х2
		MATRIX2X2_ma.M[0][0] = m[0][0]; ma.M[0][1] = m[0][1]; 
		MATRIX2X2_ma.M[1][0] = m[1][0]; ma.M[1][1] = m[1][1]; 

	 },
	 
	 Mat_Add_2X2: function(MATRIX2X2_ma, MATRIX2X2_mb, MATRIX2X2_msum)
	 {
		// Складывает две матрицы 2х2 и возвращает через msum
		MATRIX2X2_msum.M[0][0] = MATRIX2X2_ma.M[0][0]+MATRIX2X2_mb.M[0][0];
		MATRIX2X2_msum.M[0][1] = MATRIX2X2_ma.M[0][1]+MATRIX2X2_mb.M[0][1];
		MATRIX2X2_msum.M[1][0] = MATRIX2X2_ma.M[1][0]+MATRIX2X2_mb.M[1][0];
		MATRIX2X2_msum.M[1][1] = MATRIX2X2_ma.M[1][1]+MATRIX2X2_mb.M[1][1];
	 },

	 Mat_Mul_2X2: function(MATRIX2X2_ma, MATRIX2X2_mb, MATRIX2X2_mprod)
	 {
		// Перемножает две матрицы 2х2 и возвращает через mprod
		MATRIX2X2_mprod.M[0][0] = MATRIX2X2_ma.M[0][0]*MATRIX2X2_mb.M[0][0] + MATRIX2X2_ma.M[0][1]*MATRIX2X2_mb.M[1][0];
		MATRIX2X2_mprod.M[0][1] = MATRIX2X2_ma.M[0][0]*MATRIX2X2_mb.M[0][1] + MATRIX2X2_ma.M[0][1]*MATRIX2X2_mb.M[1][1];

		MATRIX2X2_mprod.M[1][0] = MATRIX2X2_ma.M[1][0]*MATRIX2X2_mb.M[0][0] + MATRIX2X2_ma.M[1][1]*MATRIX2X2_mb.M[1][0];
		MATRIX2X2_mprod.M[1][1] = MATRIX2X2_ma.M[1][0]*MATRIX2X2_mb.M[0][1] + MATRIX2X2_ma.M[1][1]*MATRIX2X2_mb.M[1][1];
	 },

	 Mat_Inverse_2X2: function(MATRIX2X2_m, MATRIX2X2_mi)
	 {
		// Вычисляет обратную матрицу 2х2 и возвращает результат через mi

		// Вычисление детериминанта
		var det = (MATRIX2X2_m.M[0][0]*m.M11 - MATRIX2X2_m.M[0][1]*MATRIX2X2_m.M[1][0]); // Mat_Det_2X2

		// Если детерминант равен 0, то неудача - возвращаем 0
		if (Math.abs(det) < Constants.Epsilon_E5)
		   return(0);
		var det_inv = 1.0/det;

		MATRIX2X2_mi.M[0][0] =  MATRIX2X2_m.M[1][1]*det_inv;
		MATRIX2X2_mi.M[0][1] = -MATRIX2X2_m.M[0][1]*det_inv;
		MATRIX2X2_mi.M[1][0] = -MATRIX2X2_m.M[1][0]*det_inv;
		MATRIX2X2_mi.M[1][1] =  MATRIX2X2_m.M[0][0]*det_inv;

		// успех
		return(1);
	 },

	 // [Tested]
	 Mat_Det_2X2: function(MATRIX2X2_m)
	 {
		// Вычисляет детерминант матрицы 2х2
		
		return(MATRIX2X2_m.M[0][0]*MATRIX2X2_m.M[1][1] - MATRIX2X2_m.M[0][1]*MATRIX2X2_m.M[1][0]);
	 },

	 //  [Tested]
	 Solve_2X2_System: function(MATRIX2X2_A, MATRIX1X2_X, MATRIX1X2_B)
	 {
		// Решает систему уравнений AX=B и вычисляет X=A(-1)*B
		// используя првило Крамера
		
		/*
		Рассмотрим систему уравнений:

		a1x + b1y = s1
		a2x + b2y = s2

		На первом шаге вычислим определитель:
		    
		    |a1 b1|
		Δ = |	  |
		    |a2 b2|

		, его называют главным определителем системы. 

		Если Δ = 0, то система имеет бесконечно много решений или несовместна (не имеет решений). 
		В этом случае правило Крамера не поможет, нужно использовать метод Гаусса.
		Если Δ ≠ 0, то система имеет единственное решение, и для нахождения корней мы должны вычислить еще два определителя:

		     |s1 b1|	    |a1 s1|
		Δx = |     | и Δy = |     |
		     |s2 b2|        |a2 s2|

		На практике вышеуказанные определители также могут обозначаться латинской буквой D.

		Корни уравнения находим по формулам:
		    Δx       Δy
		x = ---, y = ---
		     Δ        Δ
		*/
		
		// шаг 1: Вычисляем определитель А
		var det_A = this.Mat_Det_2X2(MATRIX2X2_A);

		// проверка на 0
		if (Math.abs(det_A) < Constants.Epsilon_E5)
		   return(0);

		// шаг 2: создаем матрицы-числители путем замены соответствующих столбцов матрицы А транспонированной
		// матрицей В и находим решение системы уравнений
		var work_mat = new Matrix2x2(); 

		// копируем матрицу А в рабочую
		Matrix_Copy_2x2(work_mat, MATRIX2X2_A);

		// замена столбца х
		MAT_COLUMN_SWAP_2X2(work_mat, 0, MATRIX1X2_B);

		//Вычисляем определитель новой матрицы
		var det_ABx = this.Mat_Det_2X2(work_mat);

		// решение для X00
		MATRIX1X2_X.M[0] = det_ABx/det_A;

		//////////////////////////////////////////////////////

		// копирование А в рабочую матрицу
		Matrix_Copy_2x2(work_mat, MATRIX2X2_A);

		// замена столбца у
		MAT_COLUMN_SWAP_2X2(work_mat, 1, MATRIX1X2_B);

		var det_ABy = this.Mat_Det_2X2(work_mat);

		// решение для X01
		MATRIX1X2_X.M[1] = det_ABy/det_A;

		return(1);
		
		/*
		 * var det_A = this.Mat_Det_2x2(MATRIX2X2_A);
		 * 
		 * if (Math.abs(det_A) < Constants.Epsilon_E5)
		 * 		return(0);
		 * 
		 * MATRIX1X2_X.M[0][0] = (MATRIX1X2_B.M[0][0]*MATRIX1X2_A.M[1][1] - MATRIX1X2_B.M[0][1]*MATRIX1X2_A.M[0][1])/det_A;
		 * MATRIX1X2_X.M[0][1] = (MATRIX1X2_B.M[0][1]*MATRIX1X2_A.M[0][0] - MATRIX1X2_B.M[0][0]*MATRIX1X2_A.M[1][0])/det_A;
		 * 
		 * return(1);
		 */
	 }, 

	 Mat_Add_3X3: function(MATRIX3X3_ma, MATRIX3X3_mb, MATRIX3X3_msum)
	 {
		// Складывает две матрицы 3х3 и возвращает через msum

		for (var row = 0; row < 3; row++)
    	 {
    		for (var col = 0; col < 3; col++)
        	 {
        		MATRIX3X3_msum.M[row][col] = MATRIX3X3_ma.M[row][col] + MATRIX3X3_mb.M[row][col];
        	 }

    	 }

	 },

	 Mat_Mul_VECTOR3D_3X3: function(VECTOR3D_va, MATRIX3X3_mb, VECTOR3D_mprod)
	 {
		// Перемножает вектор на матрицу 3х3 и возвращает через mprod
    	for (var col = 0; col < 3; col++)
        {
        	// Вычисляем скалярное произведение
        	var sum = 0;
        	for (var row = 0; row < 3; row++)
            {
             	sum += (VECTOR3D_va.M[row] * MATRIX3X3_mb.M[row][col]);
            }

        	// Сохраняем результат
        	VECTOR3D_mprod.M[col] = sum;
        }
	 },

	 Mat_Init_3X3: function(MATRIX3X3_ma, m00, m01, m02, m10, m11, m12, m20, m21, m22)
	 {
		// Инициализация матрицы 3х3 значениями

		MATRIX3X3_maюM00 = m00; MATRIX3X3_maюM01 = m01; MATRIX3X3_maюM02 = m02;
		MATRIX3X3_maюM10 = m10; MATRIX3X3_maюM11 = m11; MATRIX3X3_maюM12 = m12;
		MATRIX3X3_maюM20 = m20; MATRIX3X3_maюM21 = m21; MATRIX3X3_maюM22 = m22;
	 },

	 Mat_Inverse_3X3: function(MATRIX3X3_m, MATRIX3X3_mi)
	 {
		// Вычисляет обратнцю матрицу 3х3

		// сначала находим определитель
		var det = MATRIX3X3_m.M[0][0]*(MATRIX3X3_m.M[1][1]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][1]*MATRIX3X3_m.M[1][2]) - 
            	MATRIX3X3_m.M[0][1]*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[1][2]) + 
            	MATRIX3X3_m.M[0][2]*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][1] - MATRIX3X3_m.M[2][0]MATRIX3X3_m.M[1][1]);

		if (Math.abs(det) < Constants.Epsilon_E5)
   			return(0);

		var det_inv = 1.0/det;

		// 
		MATRIX3X3_mi.M[0][0] =  det_inv*(MATRIX3X3_m.M[1][1]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][1]*MATRIX3X3_m.M[1][2]);
		MATRIX3X3_mi.M[1][0] = -det_inv*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[1][2]);
		MATRIX3X3_mi.M[2][0] =  det_inv*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][1] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[1][1]);

		MATRIX3X3_mi.M[0][1] = -det_inv*(MATRIX3X3_m.M[0][1]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][1]*MATRIX3X3_m.M[0][2]);
		MATRIX3X3_mi.M[1][1] =  det_inv*(MATRIX3X3_m.M[0][0]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[0][2]);
		MATRIX3X3_mi.M[2][1] = -det_inv*(MATRIX3X3_m.M[0][0]*MATRIX3X3_m.M[2][1] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[0][1]);

		MATRIX3X3_mi.M[0][2] =  det_inv*(MATRIX3X3_m.M[0][1]*MATRIX3X3_m.M[1][2] - MATRIX3X3_m.M[1][1]*MATRIX3X3_m.M[0][2]);
		MATRIX3X3_mi.M[1][2] = -det_inv*(MATRIX3X3_m.M[0][0]*MATRIX3X3_m.M[1][2] - MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[0][2]);
		MATRIX3X3_mi.M[2][2] =  det_inv*(MATRIX3X3_m.M[0][0]*MATRIX3X3_m.M[1][1] - MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[0][1]);

		return(1);
	 },

	 Mat_Det_3X3: function(MATRIX3X3_m)
	 {
		// Вычилсяем определитель матрицы

		return(MATRIX3X3_m.M[0][0]*(MATRIX3X3_m.M[1][1]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][1]*MATRIX3X3_m.M[1][2]) - 
       		   MATRIX3X3_m.M[0][1]*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][2] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[1][2]) + 
       		   MATRIX3X3_m.M[0][2]*(MATRIX3X3_m.M[1][0]*MATRIX3X3_m.M[2][1] - MATRIX3X3_m.M[2][0]*MATRIX3X3_m.M[1][1]) );
	 },

	 Solve_3X3_System: function(MATRIX3X3_A, MATRIX1X3_X, MATRIX1X3_B)
	 {
		// Решает систему линейных уравнений А*В=Х, где А имеет размер 3х3, а матрицы В и Х размер 1х3.
		// Если решение существует, оно сохраняется в матрице Х.

		// Шаг 1: вычиляем определитель
		float det_A = Mat_Det_3X3(A);

		// проверка определителя на 0 (елси это так, то решения не существует)
		if (fabs(det_A) < EPSILON_E5)
		   return(0);

		// Шаг 2: создаем матрицы-числители путем заменыс соответствующих столбцов матрицы А транспонированной
		// матрицей В и находим решение системы
		var work_mat = new Matrix3x3();

		// решение для x /////////////////

		// копируем матрицу MATRIX3X3_A в рабочую
		MAT_COPY_3X3(work_mat, MATRIX3X3_A);

		// замена первого столбца (столбец х)
		MAT_COLUMN_SWAP_3X3(work_mat, 0, MATRIX1X3_B);

		// вычисляем определитель
		var det_ABx = this.Mat_Det_3X3(work_mat);

		// решение записываем в  MATRIX1X3_X.M[0][0]
		MATRIX1X3_X.M[0][0] = det_ABx/det_A;

		// решение для y /////////////////

		// копируем матрицу MATRIX3X3_A в рабочую
		MAT_COPY_3X3(work_mat, MATRIX3X3_A);

		// замена второго столбца (столбец y)
		MAT_COLUMN_SWAP_3X3(work_mat, 1, MATRIX1X3_B);

		// находим определитель
		var det_ABy = this.Mat_Det_3X3(work_mat);

		// решение записываем в MATRIX1X3_X.M[0][1]
		MATRIX1X3_X.M[0][1] = det_ABy/det_A;

		// решение для z /////////////////

		// копируем матрицу MATRIX3X3_A в рабочцю
		MAT_COPY_3X3(work_mat, MATRIX3X3_A);

		// заменяем третий столбец
		MAT_COLUMN_SWAP_3X3(work_mat, 2, MATRIX1X3_B);

		// вычисляем определитель
		var det_ABz = this.Mat_Det_3X3(work_mat);

		// решаем систему для MATRIX1X3_X.M[0][2]
		MATRIX1X3_X.M[0][2] = det_ABz/det_A;

		return(1);
 	},

 	Mat_Add_4X4: function(MATRIX4X4_ma, MATRIX4X4_mb, MATRIX4X4_msum)
	 {
		// Складывает две матрицы 4x4 и возвращает через msum
		for (var row = 0; row < 4; row++)
		    {
		    for (var col = 0; col < 4; col++)
		        {
		        	MATRIX4X4_msum.M[row][col] = MATRIX4X4_ma.M[row][col] + MATRIX4X4_mb.M[row][col];
		        }
		    }
	 },

	 Mat_Mul_4X4: function(MATRIX4X4_ma, MATRIX4X4_mb, MATRIX4X4_mprod)
	 {
		// Перемножает две матрицы 4х4 и возвращает через mprod

		for (var row = 0; row < 4; row++)
		{
		    for (var col = 0; col < 4; col++)
		    {
		        // вычисляем скалярное произведение строки ma 
		        // и столбца mb
		        var sum = 0;
		        for (int index = 0; index < 4; index++)
			     {
			          sum += (MATRIX4X4_ma.M[row][index] * MATRIX4X4_mb.M[index][col]);
			     }

		        MATRIX4X4_mprod.M[row][col] = sum;
		    }
		}
	 },

	 Mat_Mul_1X4_4X4: function(MATRIX1X4_ma, MATRIX4X4_mb, MATRIX1X4_mprod)
	 {
		// Перемножает матрицу размером 1х4 и 4х4 - ma*mb и возвращает через mprod

	    for (var col = 0; col < 4; col++)
	    {
		    // вычисляем скалярное произведение строки ma 
		    // и столбца mb
		    var sum = 0;
		    for (var row = 0; row < 4; row++)
		    {
		    	sum += (MATRIX1X4_ma.M[row] * MATRIX4X4_mb.M[row][col]);
		    }

		    MATRIX1X4_mprod.M[col] = sum;
	    }
	 },

	 Mat_Mul_VECTOR3D_4X4: function(VECTOR3D_va, MATRIX4X4_mb, VECTOR3D_vprod)
	 {
		// Умножает VECTOR3D на матрицу размером 
		// 4x4 matrix - ma*mb и возвращает результат через mprod
		// Функция предполагает, что вектор однородный
		// Предполагается, что w=1
		
		 for (var col = 0; col < 3; col++)
		 {
			 // вычисляет скалярное произведение строки ma и столбца mb
			 var sum = 0;

			 for (var row = 0; row < 3; row++)
			 {
			     sum += (VECTOR3D_va.M[row] * MATRIX4X4_mb.M[row][col]);
			 }

			 sum += MATRIX4X4_mb.M[row][col];    
			 VECTOR3D_vprod.M[col] = sum;
		 }
	 },

	 Mat_Mul_VECTOR3D_4X3: function(VECTOR3D_va, MATRIX4X3_mb, VECTOR3D_vprod)
	 {
		// Перемножает VECTOR3D и матрицу
		// 4x3 matrix - ma*mb и возвращает через vprod
		// Функция предполагает, что вектор однородный
		// Предполагается, что w=1

		for (var col = 0; col < 3; col++)
		{
		    // вычисляет скалярное произведение строки ma и столбца mb
		    var sum = 0;

		    for (var row = 0; row < 3; row++)
		    {
			    sum += (VECTOR3D_va.M[row] * MATRIX4X3_mb.M[row][col]);
		    }

		    sum += MATRIX4X3_mb.M[row][col];    
		    VECTOR3D_vprod.M[col] = sum;
		}
	 },

	 Mat_Mul_VECTOR4D_4X4: function(VECTOR4D_va, MATRIX4X4_mb, VECTOR4D_vprod)
	 {
		// Перемножает VECTOR4D на матрицу размером 
		// 4x4 - ma*mb и возвращает через mprod

	    for (var col = 0; col < 4; col++)
	    {
		    // вычисляем скалярное произведение строки ma и столбца mb
		    var sum = 0;

		    for (var row = 0; row < 4; row++)
		    {
			    sum += (VECTOR4D_va.M[row] * MATRIX4X4_mb.M[row][col]);
		    }

		    VECTOR4D_vprod.M[col] = sum;
	    }
	 },
	 
	 Mat_Mul_VECTOR4D_4X3: function(VECTOR4D_va, MATRIX4X4_mb, VECTOR4D_vprod)
	 {
		// Перемножает VECTOR4D и матрицу размером
		// 4x3 - ma*mb и возвращает через mprod
		// Предполагается, что последний столбец матрицы
		// mb равен [0 0 0 1] , w копируется из вектора [x y z w]

	    for (var col=0; col < 3; col++)
	    {
		    // вычисляем скалярное произведение строки ma и столбца mb
		    var sum = 0;
		    for (var row = 0; row < 4; row++)
		    {
			    sum += (VECTOR4D_va.M[row] * MATRIX4X4_mb.M[row][col]);
		    }

		    VECTOR4D_vprod.M[col] = sum;
	    }

	    VECTOR4D_vprod.M[3] = va.M[3];
	 },

	 Mat_Mul_VECTOR4D_4X3: function(VECTOR4D_va, MATRIX4X4_mb, VECTOR4D_vprod)
	 {
		// Перемножает VECTOR4D и матрицу
		// 4x3 - ma*mb и возвращает вектор через mprod
		// Предполагает, что последний столбец
		// mb равен [0 0 0 1] , w копируется из вектора [x y z w]

	    for (var col = 0; col < 3; col++)
	    {
		     // вычисляет скалярное произведение строки ma и столбца mb
		     var sum = 0;

		     for (var row = 0; row < 4; row++)
		     {
			     sum+=(VECTOR4D_va.M[row]*MATRIX4X4_mb.M[row][col]);
		     }

		     VECTOR4D_vprod.M[col] = sum;
	    }

	    VECTOR4D_vprod.M[3] = va.M[3];
	 },

	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 VECTOR3D_Print: function(va, name, outputTag)
	{
		var OutputNode = document.getElementById(outputTag);
		var val = document.createElement("p");
		var OutputString = "";
		OutputString = OutputString.concat(name, ": x=", va.x.toString(), ", y=", va.y.toString(), ", z=", va.z.toString());		
		val.innerHTML = OutputString;
		OutputNode.appendChild(val);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}


