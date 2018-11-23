/* Automatically generated code, do not edit */
/* Generated from source file: orient4d.pck */

inline int orienth_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3, const double* p4, double h0, double h1, double h2, double h3, double h4) {
    double a11;
    a11 = (p1[0] - p0[0]);
    double a12;
    a12 = (p1[1] - p0[1]);
    double a13;
    a13 = (p1[2] - p0[2]);
    double a14;
    a14 = (h1 - h0);
    double a21;
    a21 = (p2[0] - p0[0]);
    double a22;
    a22 = (p2[1] - p0[1]);
    double a23;
    a23 = (p2[2] - p0[2]);
    double a24;
    a24 = (h2 - h0);
    double a31;
    a31 = (p3[0] - p0[0]);
    double a32;
    a32 = (p3[1] - p0[1]);
    double a33;
    a33 = (p3[2] - p0[2]);
    double a34;
    a34 = (h3 - h0);
    double a41;
    a41 = (p3[0] - p0[0]);
    double a42;
    a42 = (p3[1] - p0[1]);
    double a43;
    a43 = (p3[2] - p0[2]);
    double a44;
    a44 = (h4 - h0);
    double Delta;
    Delta = ((((a11 * (((a22 * ((a33 * a44) - (a34 * a43))) - (a32 * ((a23 * a44) - (a24 * a43)))) + (a42 * ((a23 * a34) - (a24 * a33))))) - (a21 * (((a12 * ((a33 * a44) - (a34 * a43))) - (a32 * ((a13 * a44) - (a14 * a43)))) + (a42 * ((a13 * a34) - (a14 * a33)))))) + (a31 * (((a12 * ((a23 * a44) - (a24 * a43))) - (a22 * ((a13 * a44) - (a14 * a43)))) + (a42 * ((a13 * a24) - (a14 * a23)))))) - (a41 * (((a12 * ((a23 * a34) - (a24 * a33))) - (a22 * ((a13 * a34) - (a14 * a33)))) + (a32 * ((a13 * a24) - (a14 * a23))))));
    int int_tmp_result;
    double eps;
    double max1 = fabs(a11);
    if( (max1 < fabs(a21)) )
    {
        max1 = fabs(a21);
    } 
    if( (max1 < fabs(a31)) )
    {
        max1 = fabs(a31);
    } 
    if( (max1 < fabs(a41)) )
    {
        max1 = fabs(a41);
    } 
    double max2 = fabs(a12);
    if( (max2 < fabs(a22)) )
    {
        max2 = fabs(a22);
    } 
    if( (max2 < fabs(a32)) )
    {
        max2 = fabs(a32);
    } 
    if( (max2 < fabs(a42)) )
    {
        max2 = fabs(a42);
    } 
    double max3 = fabs(a13);
    if( (max3 < fabs(a14)) )
    {
        max3 = fabs(a14);
    } 
    if( (max3 < fabs(a23)) )
    {
        max3 = fabs(a23);
    } 
    if( (max3 < fabs(a24)) )
    {
        max3 = fabs(a24);
    } 
    if( (max3 < fabs(a33)) )
    {
        max3 = fabs(a33);
    } 
    if( (max3 < fabs(a34)) )
    {
        max3 = fabs(a34);
    } 
    double max4 = fabs(a23);
    if( (max4 < fabs(a24)) )
    {
        max4 = fabs(a24);
    } 
    if( (max4 < fabs(a33)) )
    {
        max4 = fabs(a33);
    } 
    if( (max4 < fabs(a34)) )
    {
        max4 = fabs(a34);
    } 
    if( (max4 < fabs(a43)) )
    {
        max4 = fabs(a43);
    } 
    if( (max4 < fabs(a44)) )
    {
        max4 = fabs(a44);
    } 
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if( (max2 < lower_bound_1) )
    {
        lower_bound_1 = max2;
    } 
    else 
    {
        if( (max2 > upper_bound_1) )
        {
            upper_bound_1 = max2;
        } 
    } 
    if( (max3 < lower_bound_1) )
    {
        lower_bound_1 = max3;
    } 
    else 
    {
        if( (max3 > upper_bound_1) )
        {
            upper_bound_1 = max3;
        } 
    } 
    if( (max4 < lower_bound_1) )
    {
        lower_bound_1 = max4;
    } 
    else 
    {
        if( (max4 > upper_bound_1) )
        {
            upper_bound_1 = max4;
        } 
    } 
    if( (lower_bound_1 < 2.89273249588395194294e-74) )
    {
        return FPG_UNCERTAIN_VALUE;
    } 
    else 
    {
        if( (upper_bound_1 > 7.23700557733225980357e+75) )
        {
            return FPG_UNCERTAIN_VALUE;
        } 
        eps = (3.17768858673611390687e-14 * (((max3 * max4) * max2) * max1));
        if( (Delta > eps) )
        {
            int_tmp_result = 1;
        } 
        else 
        {
            if( (Delta < -eps) )
            {
                int_tmp_result = -1;
            } 
            else 
            {
                return FPG_UNCERTAIN_VALUE;
            } 
        } 
    } 
    return int_tmp_result;
} 
