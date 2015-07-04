#include <cstdio>
#include <cstdlib>

// Note for all of the following delta_x = x_i_plus_1 - x_i

//****************************************************************************
// calculate_second_order_centred_finite_difference()
//****************************************************************************
double calculate_second_order_centred_finite_difference(double delta_x, double f_i_minus_1, double f_i_plus_1)
{
	double df_dx;
	df_dx = (f_i_plus_1 - f_i_minus_1) / 2.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_forth_order_centred_finite_difference()
//****************************************************************************
double calculate_forth_order_centred_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i_plus_1, double f_i_plus_2)
{
	double df_dx;
	df_dx = (f_i_minus_2 - 8.0 * f_i_minus_1 + 8.0 * f_i_plus_1 - f_i_plus_2) / 12.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_third_order_semi_centred_forward_finite_difference()
//****************************************************************************
double calculate_semi_centred_forward_finite_difference(double delta_x, double f_i_minus_1, double f_i, double f_i_plus_1, double f_i_plus_2)
{
	double df_dx;
	df_dx = (-2.0 * f_i_minus_1 - 3.0 * f_i + 6.0 * f_i_plus_1 - f_i_plus_2) / 6.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_third_order_semi_centred_backward_finite_difference()
//****************************************************************************
double calculate_semi_centred_backward_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i, double f_i_plus_1)
{
	double df_dx;
	df_dx = (f_i_minus_2 - 6.0 * f_i_minus_1 + 3.0 * f_i + 2.0 * f_i_plus_1) / 6.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_second_order_forward_finite_difference()
//****************************************************************************
double calculate_second_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2)
{
	double df_dx;
	df_dx = (-3.0 * f_i + 4.0 * f_i_plus_1 - f_i_plus_2) / 2.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_second_order_backward_finite_difference()
//****************************************************************************
double calculate_second_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i)
{
	double df_dx;
	df_dx = (f_i_minus_2 - 4.0 * f_i_minus_1 + 3.0 * f_i) / 6.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_third_order_forward_finite_difference()
//****************************************************************************
double calculate_third_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2, double f_i_plus_3)
{
	double df_dx;
	df_dx = (-11.0 * f_i + 18.0 * f_i_plus_1 - 9.0 * f_i_plus_2 + 2.0 * f_i_plus_3) / 6.0 / delta_x;
	return df_dx;
}

//****************************************************************************
// calculate_third_order_backward_finite_difference()
//****************************************************************************
double calculate_third_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i)
{
	double df_dx;
	df_dx = (-2.0 * f_i_minus_3 + 9.0 * f_i_minus_2 - 18.0 * f_i_minus_1 + 11.0 * f_i) / 6.0 / delta_x;
	return df_dx;
}
