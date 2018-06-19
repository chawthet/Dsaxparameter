package saxparam;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class TSUtils {
	
       /**
	   * Finds the minimal value in timeseries.
	   * 
	   * @param series The timeseries.
	   * @return The min value.
	   */
	  public static double min(double[] series) {
	    double min = Double.MAX_VALUE;
	    for (int i = 0; i < series.length; i++) {
	      if (min > series[i]) {
	        min = series[i];
	      }
	    }
	    return min;
	  }

	
	   /**
	   * Finds the maximal value in timeseries.
	   * 
	   * @param series The timeseries.
	   * @return The max value.
	   */
	  public static double max(double[] series) {
	    double max = Double.MIN_VALUE;
	    for (int i = 0; i < series.length; i++) {
	      if (max < series[i]) {
	        max = series[i];
	      }
	    }
	    return max;
	  }
	  /**
	   * Finds the maximal value in timeseries.
	   * 
	   * @param series The timeseries.
	   * @return The max value.
	   */
	 public static int max(int []series){
		 int max=Integer.MIN_VALUE;
		 for(int i=0;i<series.length;i++){
			 if(max<series[i]){
				 max=series[i];
			 }
		 }
		 return max;
	 }
	/**
	 * Computes the mean value of timeseries.
	 * 
	 * @param series
	 *            The timeseries.
	 * @return The mean value.
	 */
	public static double mean(double[] series) {
		double res = 0D;
		int count = 0;
		for (double tp : series) {
			res += tp;
			count += 1;

		}
		if (count > 0) {
			return res / ((Integer) count).doubleValue();
		}
		return Double.NaN;
	}

	/**
	   * Speed-optimized implementation.
	   * 
	   * @param series The timeseries.
	   * @return the standard deviation.
	   */
	  public static double stDev(double[] series) {
	    double num0 = 0D;
	    double sum = 0D;
	    int count = 0;
	    for (double tp : series) {
	      num0 = num0 + tp * tp;
	      sum = sum + tp;
	      count += 1;
	    }
	    double len = ((Integer) count).doubleValue();
	    return Math.sqrt((len * num0 - sum * sum) / (len * (len - 1)));
	  }

	 /*
	  * Normalize the data between 0 and 1
	  * z_i= x_i - min(x)/ max(x)-min(x) 
	  */
	  public static double[]scaleNorm(double[] series){
		  double max= max(series);
		  double min=min(series);
		  double[] result=new double[series.length];
		  for(int i=0;i< series.length;i++){
			  result[i]=(series[i]-min)/(max-min);
		  }
		  return result;
	  }
	  
	/**
	 * Speed-optimized Z-Normalize routine, doesn't care about normalization
	 * threshold.
	 * 
	 * @param series
	 *            The timeseries.
	 * @param normalizationThreshold
	 *            the zNormalization threshold value.
	 * @return Z-normalized time-series.
	 */
	public static double[] znorm(double[] series, double normalizationThreshold) {
		double[] res = new double[series.length];
		double mean = mean(series);
		double sd = stDev(series);
		if (sd < normalizationThreshold) {
			return series.clone();
		}
		for (int i = 0; i < res.length; i++) {
			res[i] = (series[i] - mean) / sd;
		}
		return res;
	}
	
	public static ArrayList<Double>znorm(ArrayList<Double>series, double normalizationThreshold){
		ArrayList<Double>res=new ArrayList<Double>();
		double []tmpSeries=new double[series.size()];
		for(int i=0;i<tmpSeries.length;i++){
			tmpSeries[i]=series.get(i);
		}
		double mean=mean(tmpSeries);
		double sd=stDev(tmpSeries);
		if(sd<normalizationThreshold){
			return series;
		}
		for(int i=0;i<tmpSeries.length;i++){
			res.add((tmpSeries[i]-mean)/sd);
		}
		return res;
	}
	public static double[] znorm(double[] series) {
		double[] res = new double[series.length];
		double mean = mean(series);
		double sd = stDev(series);
		
		for (int i = 0; i < res.length; i++) {
			res[i] = (series[i] - mean) / sd;
		}
		return res;
	}
	
	/*public static double[] znorm(double[] series) {
		double[] res = new double[series.length];
		double mean = mean(series);
		double sd = stDev(series);		
		for (int i = 0; i < res.length; i++) {
			res[i] = (series[i] - mean) / sd;
		}
		return res;
	}
	 */
	
	// DataLoad
	public static List<Double> dataLoad(String filename) throws Exception {
		Charset charset = StandardCharsets.UTF_8;
		Path file = Paths.get(filename);
		if (!(Files.exists(file))) {
			throw new Exception("unable to load data - data source not found.");
		}
		List<Double> dList = new ArrayList<Double>();
		BufferedReader reader;
		try {
			reader = Files.newBufferedReader(file, charset);
			String line = null;
			while ((line = reader.readLine()) != null) {
				dList.add(Double.parseDouble(line));
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return dList;
	}
	
	
	//check the arraysize and compare PowerOfTwo
	public static int checkPowerOfTwo(int TSLength){
		int []powerOfTwo={16,32,64,128,256,512,1024,2048};
		int diffLen=0;
		for(int i=0;i<powerOfTwo.length;i++){
			if(TSLength <= powerOfTwo[i]){
				diffLen=powerOfTwo[i]-TSLength;
				break;
			}
		}
		return diffLen;
	}
	//Filling Array with Zero for creating PowerOfTwo arraySize
	public static double[] fillZeroarray(int diffLen, double[]orgArray){
		double[]newArray=new double[diffLen+orgArray.length];
		double[]zeroArr=new double[diffLen];
		Arrays.fill(zeroArr, 0);
		System.arraycopy(orgArray, 0, newArray, 0, orgArray.length);
		System.arraycopy(zeroArr, 0, newArray, orgArray.length, zeroArr.length);
		return newArray;
	}
	
	
	//Finding local Minima I
	public static int findLocalMinima(double []arr, int lowIndex, int highIndex){
		if(lowIndex > highIndex)
			return -1;
		if(lowIndex == highIndex)
			return lowIndex;
		if(lowIndex+1 == highIndex){
			if (arr[lowIndex] < arr[highIndex]){
				return lowIndex;
			}
			return highIndex;
		}
		int midIndex=(lowIndex + highIndex)/2;
		if(arr[midIndex]<arr[midIndex-1] && arr[midIndex]<= arr[midIndex+1])
			return midIndex;
		if(arr[midIndex]>arr[midIndex+1])
			return findLocalMinima(arr,midIndex,highIndex);
		else
			return findLocalMinima(arr,lowIndex,midIndex);
	}
	
	
	public static class Index{
		private double value;
		private int index;
		public Index(){
			
		}
		public Index(double value, int index){
			this.value=value;
			this.index=index;			
		}
		public double getValue() {
			return value;
		}
		public void setValue(double value) {
			this.value = value;
		}
		public int getIndex() {
			return index;
		}
		public void setIndex(int index) {
			this.index = index;
		}
		
	}
	public static class IndexComparator implements Comparator<Object>{

		@Override
		public int compare(Object o1, Object o2) {
			Index i1=new Index();
			i1=(Index) o1;
			Index i2=new Index();
			i2=(Index) o2;
			return (int) Math.abs(i2.getValue() - Math.abs(i1.getValue()));			
		}
		
	}
	
	
	//Finding local Minima II
	public static int findLocalMinima(double []arr){
	
		
		Index[]arr_Idx=new Index[arr.length];
		for(int i=0;i< arr.length;i++){
			arr_Idx[i]=new Index(arr[i],i);
		}
		
		Index[]aux_Idx=new Index[arr.length];
		
		
		Arrays.sort(arr_Idx,new IndexComparator());
		
		double maxVal=arr_Idx[arr.length-1].getValue();
		
		for(int i=0;i< arr.length;i++){
			double temp=maxVal-arr[i];
			aux_Idx[i]=new Index(temp,i);
			
		}
		Arrays.sort(aux_Idx,new IndexComparator());
		return aux_Idx[arr.length-1].getIndex();		
	}			
	
	//finding local minima
		public static int lminima(double[]tArray){
			for(int i=1;i< tArray.length-1;i++){
				if(tArray[i-1]>tArray[i] && tArray[i+1]>tArray[i]){
					return i;
				}
				else continue;			
			}
		return 0;
		}
}
