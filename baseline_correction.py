import numpy as np
import matplotlib.pyplot as plt

# def find_baseline_end_region(x, y, threshold=0.1, min_points=3, method="mean"):
#     """
#     Finds the flat baseline region at the end of the dataset.
    
#     Parameters:
#         x (array): The x-values of the data.
#         y (array): The y-values of the data.
#         threshold (float): Maximum difference between consecutive y-values to consider flat.
#         min_points (int): Minimum number of points required in the flat region.
#         method (str): How to calculate the baseline value ("mean" or "median").
    
#     Returns:
#         float: The baseline value for the flat region.
#     """
#     # Calculate the difference between consecutive y-values
#     dy = np.abs(np.diff(y))
    
#     # Find indices where the difference is below the threshold
#     flat_indices = np.where(dy < threshold)[0]
    
#     # Check for consecutive indices forming a flat region near the end
#     end_flat_indices = []
#     for i in range(len(flat_indices) - min_points + 1):
#         if np.all(np.diff(flat_indices[i:i + min_points]) == 1):
#             end_flat_indices = flat_indices[i:i + min_points]
#             break

#     if len(end_flat_indices) == 0:
#         raise ValueError("No sufficient flat region found at the end. Try adjusting 'threshold' or 'min_points'.")
    
#     # Extract the y-values for the flat region
#     baseline_region_y = y[end_flat_indices[0]:end_flat_indices[-1] + 2]

#     # Calculate the baseline value using the specified method
#     if method == "mean":
#         baseline = np.mean(baseline_region_y)
#     elif method == "median":
#         baseline = np.median(baseline_region_y)
#     else:
#         raise ValueError("Invalid method. Use 'mean' or 'median'.")
    
#     return baseline, end_flat_indices

# def find_baseline_end_region(x, y, threshold=0.05, min_points=3, method="mean", end_region_width=5):
#     """
#     Identify the baseline as the flat region near the end of the data.
#     """
#     x = np.array(x)
#     y = np.array(y)
#     dy = np.abs(np.diff(y))

#     print("Differences between y values (dy):", dy)

#     flat_indices = np.where(dy < threshold)[0]
#     print("Flat indices (where dy < threshold):", flat_indices)

#     # end_flat_indices = [idx for idx in flat_indices if idx > len(y) - min_points - 1]
#     # Relaxed condition for "end flat indices"
#     end_flat_indices = [idx for idx in flat_indices if idx >= len(y) - end_region_width]
#     print("End flat indices (near end of dataset):", end_flat_indices)

#     if len(end_flat_indices) < min_points:
#         raise ValueError("No sufficient flat region found at the end. Try adjusting 'threshold' or 'min_points'.")

#     flat_y = y[end_flat_indices]
#     if method == "mean":
#         baseline = np.mean(flat_y)
#     elif method == "median":
#         baseline = np.median(flat_y)
#     else:
#         raise ValueError("Invalid method. Use 'mean' or 'median'.")
    
#     return baseline


# # Test the function
# y = np.array([6.39, 7.19, 8.17, 7.52, 8.01, 8.74, 8.18, 8.86, 9.02, 9.98,
#      9.4, 10.48, 9.59, 11.22, 11.40, 13.22, 14.68, 14.12,
#      15.68, 16.75, 17.508, 19.45, 20.06, 24.002, 25.969,
#      28.237, 23.788, 4.534, 2.417, 2.734])
# y = np.array([
#     20.64115698, 22.36628989, 22.64895175, 23.94208498, 23.5669892, 
#     24.65534819, 25.3661616, 25.1744941, 26.52510513, 28.35869341, 
#     10.95475229, 6.17815451, 6.27072533, 6.11349406, 5.79806452, 
#     5.28778689, 5.22739726, 5.63908304, 5.47278574, 5.42443154, 
#     5.57195391, 5.43084224, 5.40802596, 5.34045332, 5.08130624, 
#     5.30625598, 5.38383312, 5.34101267, 5.74487288, 4.85098167
# ])
# x = np.arange(len(y))

# # Detect the flat baseline region at the end
# # baseline = find_baseline_end_region(x, y, threshold=2, min_points=3, method="mean", end_region_width=5)
# baseline, flat_region_indices = find_baseline_end_region(x, y, threshold=1, min_points=3, method="mean")
# print(f"Baseline: {baseline}")
# print(f"Flat region indices: {flat_region_indices}")


# import matplotlib.pyplot as plt
# # Plot the data
# plt.plot(x, y, label="Data")
# # Highlight the flat region (baseline)
# flat_indices = np.arange(len(y) - 3, len(y))
# plt.axhline(y=baseline, color="red", linestyle="--", label="Baseline")
# # Mark flat points
# plt.scatter([x[i] for i in flat_indices], [y[i] for i in flat_indices], color="green", label="Flat Region")
# plt.legend()
# plt.show()



def moving_average(y, window_size):
    """
    Compute the moving average of the input data y using a specified window size.
    
    Parameters:
    y (numpy array): The data for which the moving average is to be calculated.
    window_size (int): The size of the moving window.
    
    Returns:
    numpy array: The moving average of the input data.
    """
    return np.convolve(y, np.ones(window_size)/window_size, mode='valid')



# def find_flat_region_at_end(y, threshold=0.1, window_size=5, tail_length=10):
#     """
#     Find a flat region at the end of the data where the values are within a certain threshold.
    
#     Parameters:
#     y (numpy array): The input data to find the flat region.
#     threshold (float): The allowable variation in the data to consider it "flat."
#     window_size (int): The size of the moving average window.
#     tail_length (int): The number of points from the end to analyze for flatness.
    
#     Returns:
#     numpy array: The indices of the flat region toward the end.
#     """
#     # Focus on the last `tail_length` values
#     y_tail = y[-tail_length:]
    
#     # Compute the moving average for the tail part of the data
#     smoothed_y_tail = moving_average(y_tail, window_size)
    
#     # Compute the differences between adjacent values in the smoothed tail data
#     diff = np.abs(np.diff(smoothed_y_tail))
    
#     # Find regions where the difference is smaller than the threshold
#     flat_indices = np.where(diff < threshold)[0]
    
#     # Convert indices back to the original indices
#     flat_indices_original = len(y) - tail_length + flat_indices
    
#     return flat_indices_original

# def find_baseline(y, threshold=0.1, window_size=5, tail_length=10):
#     flat_indices_y = find_flat_region_at_end(y, threshold, window_size, tail_length)
#     flat_points_y = y[flat_indices_y]
#     # print("Flat points y:", flat_points_y)
#     min_flat_y = np.min(flat_points_y)
#     mean_flat_y = np.mean(flat_points_y)
#     # print("Min flat y:", min_flat_y)
#     # print("Mean flat y:", mean_flat_y)
#     diff = np.abs(mean_flat_y - min_flat_y)
#     # print("Difference:", diff)
#     if diff < 2:
#         baseline_aray = flat_points_y
#         baseline_indices = flat_indices_y
#     elif diff >= 2:
#         baseline_aray = np.array([min_flat_y])
#         baseline_indices = np.where(y == min_flat_y)[0]
#     return baseline_aray, baseline_indices



# Test the function
y1 = np.array([6.39, 7.19, 8.17, 7.52, 8.01, 8.74, 8.18, 8.86, 9.02, 9.98,
     9.4, 10.48, 9.59, 11.22, 11.40, 13.22, 14.68, 14.12,
     15.68, 16.75, 17.508, 19.45, 20.06, 24.002, 25.969,
     28.237, 23.788, 4.534, 2.417, 2.734])
y2 = np.array([
    20.64115698, 22.36628989, 22.64895175, 23.94208498, 23.5669892, 
    24.65534819, 25.3661616, 25.1744941, 26.52510513, 28.35869341, 
    10.95475229, 6.17815451, 6.27072533, 6.11349406, 5.79806452, 
    5.28778689, 5.22739726, 5.63908304, 5.47278574, 5.42443154, 
    5.57195391, 5.43084224, 5.40802596, 5.34045332, 5.08130624, 
    5.30625598, 5.38383312, 5.34101267, 5.74487288, 4.85098167
])

x = np.arange(len(y1))

y3 = np.array([
    19.30883442, 20.07757335, 18.96119217, 20.86382886, 20.48170201, 21.36438848,
    21.63018098, 21.78063575, 23.27558741, 24.93535953, 23.59854513, 26.10346948,
    25.56415321, 28.99915235, 32.45695705, 36.33465568, 41.47629481, 42.45185953,
    33.99725735, 8.00382278, 6.39998534, 6.1643036, 6.30447165, 6.42192978,
    6.21152115, 6.16177767, 6.5173214, 6.06197459, 6.13204433, 5.93316579
])

y4 = np.array([
    22.1444604, 23.01215988, 23.01305708, 26.04791137, 26.28603335, 28.16893052,
    27.90362039, 29.36587501, 31.1409771, 35.46440521, 32.18174645, 37.97447576,
    35.16580886, 27.20743578, 5.55349661, 5.51783765, 5.64663869, 5.40597492,
    5.4658837, 5.57264111, 5.34404792, 5.66103638, 5.33349183, 5.52470221,
    5.62855117, 5.68862557, 5.49138078, 5.69161193, 5.60824, 5.16557523
])

y5 = np.array([
    24.49296548, 25.46627653, 27.47677258, 28.44274262, 29.84240352, 32.33579699,
    34.76878934, 35.68491466, 40.22614988, 40.30516193, 21.12344565, 5.58817917,
    6.13746806, 5.9098371, 5.85587476, 4.93791832, 5.23936212, 5.55480845,
    5.38742329, 5.42908318, 5.44059027, 5.26325839, 5.45694016, 5.32882292,
    5.38560345, 5.04186603, 5.09603059, 5.32197719, 5.20021575, 4.67376361
])

y6 = np.array([
    6.86746039, 6.86174839, 7.22660142, 7.18102129, 6.69747302, 7.40070116,
    8.28359464, 6.86598553, 8.83833186, 8.46581126, 8.97176758, 9.1037164,
    8.74981295, 9.48991126, 9.94919109, 11.67198452, 12.0182097, 12.15037258,
    13.99077146, 13.37034899, 15.32036894, 15.31188849, 16.1815145, 18.00098574,
    20.59605422, 22.6842125, 20.00806372, 6.1574421, 1.24400165, 1.30416514
])   # event Nbr 15736 of 1005 file

y7 = np.array([
    7.99115087, 8.8731933, 8.62346948, 9.3155056, 9.10757759, 9.82433461,
    9.82761015, 10.21501279, 10.865725, 10.83596941, 11.11509707, 11.70866108,
    11.44794399, 12.60570605, 13.85539448, 14.68005782, 16.15993464, 15.47958455,
    17.659119, 18.28436637, 19.52873129, 20.00990523, 22.19400068, 22.51017344,
    27.32119477, 28.96951372, 31.4497906, 32.6843946, 21.50476169, 1.18248749
])   # event Nbr 415653 of 1005 file

def fn_baseline_correction(y, window_size=2, threshold=0.2):
    dy = np.abs(np.diff(y))
    print("Differences between y values (dy):", dy)
    max_diff = np.max(dy)
    print("Max difference:", max_diff)
    max_diff_index = np.where(dy == max_diff)[0]
    print("Max difference index:", max_diff_index)
    print("y tail section: ", y[max_diff_index[0]:])
    y_tail = y[max_diff_index[0]:]
    # Compute the moving average for the tail part of the data
    smoothed_y_tail = moving_average(y_tail, window_size)
    print("Smoothed y tail:", smoothed_y_tail)
    # Compute the differences between adjacent values in the smoothed tail data
    diff = np.abs(np.diff(smoothed_y_tail))
    print("Differences in smoothed y tail:", diff)
    if max_diff_index[0] > 25:
        threshold = 5.0
    else:
        threshold = threshold
    # Find regions where the difference is smaller than the threshold
    flat_indices = np.where(diff < threshold)[0]
    print("Flat indices shape: ", flat_indices.shape)
    # for safety check, if flat_indices is empty, set it to 0
    if flat_indices.size == 0:
        flat_indices = np.array([0])
    else:
        pass 
    print("Flat indices:", flat_indices)
    baseline_corrected_y = y
    baseline_corrected_y[max_diff_index[0]+flat_indices[0]:]= 0.0
    print("Baseline corrected y:", baseline_corrected_y)
    return baseline_corrected_y


# print("Baseline corrected y3:", fn_baseline_correction(y3))

# # Set the window size and threshold parameters
window_size = 2
threshold = 0.3  # You can adjust this to define what you consider "flat"

# Plotting the results
plt.figure(figsize=(20, 15))

plt.subplot(3, 3, 1)
plt.plot(x, y1, label='y1')
flat_indices_y1 = np.where(fn_baseline_correction(y1, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y1], y1[flat_indices_y1], 'ro', label='Flat Region')
plt.legend()
plt.title('y1 with Flat Region at the End')

plt.subplot(3, 3, 2)
plt.plot(x, y2, label='y2')
flat_indices_y2 = np.where(fn_baseline_correction(y2, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y2], y2[flat_indices_y2], 'ro', label='Flat Region')
plt.legend()
plt.title('y2 with Flat Region at the End')

plt.subplot(3, 3, 3)
plt.plot(x, y3, label='y3')
flat_indices_y3 = np.where(fn_baseline_correction(y3, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y3], y3[flat_indices_y3], 'ro', label='Flat Region')
plt.legend()
plt.title('y3 with Flat Region at the End')

plt.subplot(3, 3, 4)
plt.plot(x, y4, label='y4')
flat_indices_y4 = np.where(fn_baseline_correction(y4, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y4], y4[flat_indices_y4], 'ro', label='Flat Region')
plt.legend()
plt.title('y4 with Flat Region at the End')

plt.subplot(3, 3, 5)
plt.plot(x, y5, label='y5')
flat_indices_y5 = np.where(fn_baseline_correction(y5, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y5], y5[flat_indices_y5], 'ro', label='Flat Region')
plt.legend()
plt.title('y5 with Flat Region at the End')

plt.subplot(3, 3, 6)
plt.plot(x, y6, label='y6')
flat_indices_y6 = np.where(fn_baseline_correction(y6, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y6], y6[flat_indices_y6], 'ro', label='Flat Region')
plt.legend()
plt.title('y6 with Flat Region at the End')

plt.subplot(3, 3, 7)
plt.plot(x, y7, label='y7')
flat_indices_y7 = np.where(fn_baseline_correction(y7, window_size, threshold)==0)[0]
plt.plot(x[flat_indices_y7], y7[flat_indices_y7], 'ro', label='Flat Region')
plt.legend()
plt.title('y7 with Flat Region at the End')


plt.tight_layout()
plt.show()


