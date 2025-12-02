import numpy as np
import time
import matplotlib.pyplot as plt

TABLE_SIZE = 128 


taylor_coeff = [ 
  1.0000000000000000,
  0.33333333333333331,
  0.13333333333333333,
  0.053968253968253971,
  0.021869488536155203,
  0.0088632355299021973,
  0.0035921280365724811,
  0.0014558343870513183,
  0.00059002744094558595,
  0.00023912911424355248,
  9.6915379569294509e-05,
  3.9278323883316833e-05,
  1.5918905069328964e-05,
  6.4516892156554306e-06
]

def taylor_tan(x):
    res = taylor_coeff[-1]
    x2 = x**2

    for c in reversed(taylor_coeff[:-1]):
        res = x2 * res + c
    
    return res * x


def get_addon(x):
    res = []

    for item in x:
        y_wanted = np.tan(item)

        lower = item 
        upper = np.pi
        midpoint = (upper-lower) /2
        eval_y = -1000

        while abs(y_wanted - eval_y) > 10e-16:
            eval_y = taylor_tan(midpoint)

            if y_wanted - eval_y > 0:
                lower = midpoint
            else:
                upper = midpoint 

            midpoint = (upper-lower) / 2 + lower
        
        res.append(midpoint)

    return np.array(res) - x


def get_index(x):
    i = (int) (1 / (np.pi/4+1/TABLE_SIZE-x))-1
    return i


def get_table(func, index_function, start_range, end_range):
    res = []
    THRESHOLD = 10e-16

    lower = start_range
    upper = end_range
    midpoint = (upper - lower) / 2 + lower 

    start_section = lower 
    end_section = upper 

    cur_index = 0


    while cur_index < TABLE_SIZE-1:
        start_section = lower

        while upper - lower > THRESHOLD:
            if index_function(midpoint) == index_function(lower):
                lower = midpoint
            else:
                upper = midpoint

            midpoint = (upper - lower) / 2 + lower 


        #print(index_function(lower), index_function(upper), cur_index)
        #time.sleep(0.01)


        end_section = midpoint

        dist = (end_section - start_section) / 2
        mid_section = dist + start_section 

        start_value = func(np.array([start_section]))[0]
        mid_value = func(np.array([mid_section]))[0]

        intersec = mid_value
        lower = end_section 
        upper = end_range

        if get_index(mid_section) == cur_index:
            res.append(intersec)
            cur_index += 1


    return res 


def apply_table(x, table):
    i = get_index(x)
    # print(table[i][1])
    return table[min(i, len(table)-1)]


def corrected_taylor(x_vals):
    x_vals = x_vals / 2 
    y_vals = taylor_tan(x_vals)
    return 2 * y_vals / (1 - y_vals **  2)


x_vals = np.linspace(0, np.pi/4, 10000)
indexes = np.array([get_index(x) for x in x_vals])
addon_vec = get_addon(x_vals)
adjustment_table = get_table(get_addon, get_index, 0, np.pi/4)


print("Max of Index: ", np.max(indexes))

"""
print("double intersec_lookup[256] = {") #}

for i, v in enumerate(adjustment_table):
    print(f"{v}", end='')
    if i != 255:
        print(",")
    
print("}")



print("double slope_lookup[256] = {") #}

for i, v in enumerate(adjustment_table):
    print(f"{v}", end='')
    if i != 255:
        print(",")
    
print("}")
"""
# y_true = np.tan(x_vals)

plt.plot(x_vals, corrected_taylor(x_vals) - np.tan(x_vals))

# y_tayl = taylor_tan(x_vals)
# plt.plot(x_vals, corrected_taylor(x_vals), label="Corrected")
# plt.plot(x_vals, np.tan(x_vals), label="True")

plt.legend()
plt.show()

