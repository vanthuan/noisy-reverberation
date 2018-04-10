function [z] = hertz2bark(h),

z = 13* atan(0.00076 * h) + 3.5 * atan((h/7500).^2);
end