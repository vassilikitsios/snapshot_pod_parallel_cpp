double calculate_second_order_centred_finite_difference(double delta_x, double f_i_minus_1, double f_i_plus_1);
double calculate_forth_order_centred_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i_plus_1, double f_i_plus_2);
double calculate_semi_centred_forward_finite_difference(double delta_x, double f_i_minus_1, double f_i, double f_i_plus_1, double f_i_plus_2);
double calculate_semi_centred_backward_finite_difference(double delta_x, double f_i_minus_2, double f_i_minus_1, double f_i, double f_i_plus_1);
double calculate_second_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2);
double calculate_second_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i);
double calculate_third_order_forward_finite_difference(double delta_x, double f_i, double f_i_plus_1, double f_i_plus_2, double f_i_plus_3);
double calculate_third_order_backward_finite_difference(double delta_x, double f_i_minus_3, double f_i_minus_2, double f_i_minus_1, double f_i);
