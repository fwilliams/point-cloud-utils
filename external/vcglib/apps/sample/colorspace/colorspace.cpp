
// Standard headers
#include <iostream>

// VCG headers
#include <vcg/space/color4.h>
#include <vcg/space/colorspace.h>

using namespace std;

typedef vcg::ColorSpace<double> ColorSpace;

int main(int argc,char ** argv)
{
	// Some Color Space processing examples

	cout << endl;

	if (argc != 4)
	{
		cout << "  Usage: colorspace <r> <g> <b>" << endl << endl;
		cout << "         <r> <g> <b> : RGB triplet (0-255)" << endl << endl;
		cout << "         Note: The RGB triplet is assumed in the sRGB space (D65 illuminant).";
		cout << endl << endl;
		exit(-2);
	}

	double r = static_cast<double>(atof(argv[1]));
	double g = static_cast<double>(atof(argv[2]));
	double b = static_cast<double>(atof(argv[3]));

	// RGB components have to be in the range [0.0, 1.0]
	vcg::Color4<double> color(r/255.0, g/255.0, b/255.0, 0.0);


	// RGB --> RGB (RGB space changing)
	///////////////////////////////////////////////

	cout << "  * RGB --> RGB conversion" << endl << endl;

	vcg::Color4<double> rgb = ColorSpace::RGBtoRGB(color, ColorSpace::SRGB, ColorSpace::PAL_RGB);
	rgb *= 255.0;
	cout << "  RGB (PAL/SECAM): " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;
	
	rgb = ColorSpace::RGBtoRGB(color, ColorSpace::SRGB, ColorSpace::WIDE_GAMUT);
	rgb *= 255.0;
	cout << "  RGB (Wide Gamut): " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl << endl;

	
	// RGB <--> HSL (Hue, Saturation, Lightness)
	///////////////////////////////////////////////

	cout << "  * RGB <--> HSL conversion" << endl << endl;

	vcg::Color4<double> hsl = ColorSpace::RGBtoHSL(color);
	cout << "  RGB --> HSL: " << hsl[0]*360.0 << "° " << hsl[1]*100.0 << "% " << hsl[2]*100.0 << "% " << endl;

	rgb = ColorSpace::HSLtoRGB(hsl);
	rgb *= 255.0;
	cout << "  HSL --> RGB: " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl << endl;


	// RGB <--> HSV (Hue, Saturation, Value)
	///////////////////////////////////////////////

	cout << "  * RGB <--> HSV conversion" << endl << endl;

	vcg::Color4<double> hsv = ColorSpace::RGBtoHSV(color);
	cout << "  RGB --> HSV: " << hsv[0]*360.0 << "° " << hsv[1]*100.0 << "% " << hsv[2]*100.0 << "% " <<  endl;

	rgb = ColorSpace::HSVtoRGB(hsv);
	rgb *= 255.0;
	cout << "  HSV --> RGB: " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl << endl;


	// RGB --> CIELab
	///////////////////////////////////////////////

	cout << "  * RGB <--> Lab conversion" << endl << endl;

	vcg::Color4<double> xyzD65 = ColorSpace::RGBtoXYZ(color, ColorSpace::SRGB, ColorSpace::VCG_ILLUMINANT_D65); 
	vcg::Color4<double> lab = ColorSpace::XYZtoCIELab(xyzD65, ColorSpace::refIlluminant(ColorSpace::SRGB));
	cout << "  RGB --> CIELab: " << lab[0] << " " << lab[1] << " " << lab[2] << endl;

	vcg::Color4<double> xyz = ColorSpace::CIELabtoXYZ(lab, ColorSpace::refIlluminant(ColorSpace::SRGB));
	rgb = ColorSpace::XYZtoRGB(xyz, ColorSpace::refIlluminant(ColorSpace::SRGB), ColorSpace::SRGB);
	rgb *= 255.0;
	cout << "  CIELab --> RGB: " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl << endl;



	// RGB <--> XYZ
	///////////////////////////////////////////////

	cout << "  * RGB <--> XYZ conversion" << endl << endl;

	// RGB --> XYZ (D65)
	cout << "  RGB --> XYZ (D65): " << xyzD65[0] << " " << xyzD65[1] << " " << xyzD65[2] << endl;

	// XYZ (D65) --> RGB
	rgb = ColorSpace::XYZtoRGB(xyzD65, ColorSpace::VCG_ILLUMINANT_D65, ColorSpace::SRGB);
	rgb *= 255.0;
	cout << "  XYZ (D65) --> RGB: " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;

	// RGB --> XYZ (D50)
	vcg::Color4<double> xyzD50 = ColorSpace::RGBtoXYZ(color, ColorSpace::SRGB, ColorSpace::VCG_ILLUMINANT_D50);
	cout << "  RGB --> XYZ (D50): " << xyzD50[0] << " " << xyzD50[1] << " " << xyzD50[2] << endl;

	// XYZ (D50) --> RGB
	rgb = ColorSpace::XYZtoRGB(xyzD50, ColorSpace::VCG_ILLUMINANT_D50, ColorSpace::SRGB);
	rgb *= 255.0;
	cout << "  XYZ (D50) --> RGB: " << rgb[0] << " " << rgb[1] << " " << rgb[2] << endl;

	// Direct way
	xyz = ColorSpace::chromaticAdaptation(xyzD65, ColorSpace::VCG_ILLUMINANT_D65, ColorSpace::VCG_ILLUMINANT_D50);
	cout << "  XYZ (D65 --> D50): " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl << endl;

	return 0;
}

