# Statistics

mean(x) = sum(x)/length(x)

# Streaming statistics, textbook formulas
update_mean(x_bar, x_new, n) = x_bar + (x_new - x_bar)/n
update_var(s2, x_bar, x_bar_old, x_new, n) =
    s2 + ((x_new - x_bar)*(x_new - x_bar_old) - s2)/n

# Streaming statistics, Welford's online algorithm
update_var_welford(S, x_bar, x_bar_old, x_new) = 
    S + (x_new - x_bar_old)*(x_new - x_bar)

# Streaming weighted statistics, textbook formulas
update_wmean(x_bar, w, x_new, w_new) =
    x_bar + (x_new - x_bar)*(w_new/(w + w_new))
update_wvar(s2, x_bar, w, x_new, x_bar_new, w_new) =
    (s2*w + w_new*(x_new - x_bar)*(x_new - x_bar_new))/(w + w_new)
