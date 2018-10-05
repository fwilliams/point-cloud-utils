#ifndef __VCGLIB_PRECISION
#define __VCGLIB_PRECISION

namespace vcg
{
	namespace tri
	{
		namespace io
		{
			template<typename SCALAR>
			struct Precision
			{
				static int digits() {return 0;}
                static const char* typeName() {return "";}
			};
			
			template<>
			struct Precision<float>
			{
				static int digits() {return 7;}
                static const char* typeName() {return "float";}

			};
			
			template<>
			struct Precision<double>
			{
				static int digits() {return 16;}
                static const char* typeName() {return "double";}
			};
		}
	}
}


#endif