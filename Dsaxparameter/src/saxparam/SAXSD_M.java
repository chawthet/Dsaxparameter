package saxparam;

import java.io.IOException;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.Frequency;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.Precision;

import net.seninp.jmotif.sax.SAXException;
import net.seninp.jmotif.sax.TSProcessor;
import net.seninp.jmotif.sax.alphabet.Alphabet;
import net.seninp.jmotif.sax.alphabet.NormalAlphabet;
import saxparam.TSUtils;

//SAXSD_M
public class SAXSD_M {
	
	// simple class to model instances (class + features)
		static class sampleSeries {
			List<Double> Attributes = new ArrayList<Double>();
			int cName;

			public sampleSeries(int cName, ArrayList<Double> Attribute) {
				this.cName = cName;
				this.Attributes = Attribute;
			}
		}	
	public static List<sampleSeries> dataLoad(String filename) {
		Path file = Paths.get(filename);
		// java 8: Stream class
		Stream<String> lines;
		int cname = 0;
		List<sampleSeries> sSeriesList = new ArrayList<sampleSeries>();
		try {
			lines = Files.lines(file, StandardCharsets.UTF_8);
			for (String line : (Iterable<String>) lines::iterator) {
				StringTokenizer stk = new StringTokenizer(line, ",");
				cname = Integer.parseInt(stk.nextToken());
				ArrayList<Double> sList = new ArrayList<Double>();
				while (stk.hasMoreTokens()) {
					sList.add(Double.parseDouble(stk.nextToken()));
				}
				sSeriesList.add(new sampleSeries(cname, sList));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sSeriesList;
	}	
	//Hierarchical Mean segmentation
		public static ArrayList<Integer> hmean_seg(double[] ts, int sIdx) {
			double mean = TSUtils.mean(ts);
			// System.out.println(mean);
			ArrayList<Integer> meanLessList = new ArrayList<>();
			ArrayList<Integer> meanGreatList = new ArrayList<>();
			for (int i = 0; i < ts.length; i++) {
				if (ts[i] < mean) {
					int tmp = i;
					meanLessList.add(tmp+sIdx);
				} else {
					int tmp = i;
					meanGreatList.add(tmp+sIdx);
				}
			}
			// display
			/*
			 * for(int i=0;i< meanLessList.size();i++){
			 * System.out.print(meanLessList.get(i)+" , "); } System.out.println();
			 * 
			 * for(int i=0;i< meanGreatList.size();i++){
			 * System.out.print(meanGreatList.get(i)+" , "); } System.out.println();
			 */
			// ArrayList<Integer>RmeanLessList=new ArrayList<>();
			Set<Integer> RmeanLessList = new HashSet<Integer>();
			Set<Integer> RmeanGreatList = new HashSet<Integer>();
			// ArrayList<Integer>RmeanGreatList=new ArrayList<>();
			RmeanLessList.add(meanLessList.get(0));
			for (int i = 1; i < meanLessList.size() - 1; i++) {
				if (meanLessList.get(i + 1) - meanLessList.get(i) > 3) {
					RmeanLessList.add(meanLessList.get(i));
					RmeanLessList.add(meanLessList.get(i + 1));
				} else
					continue;
			}
			// if(RmeanLessList.get(RmeanLessList.size()-1)!=
			// meanLessList.get(meanLessList.size()-1))
			// {
			RmeanLessList.add(meanLessList.get(meanLessList.size() - 1));
			// }
			RmeanGreatList.add(meanGreatList.get(0));
			for (int i = 1; i < meanGreatList.size() - 1; i++) {
				if (meanGreatList.get(i + 1) - meanGreatList.get(i) > 3) {
					RmeanGreatList.add(meanGreatList.get(i));
					RmeanGreatList.add(meanGreatList.get(i + 1));
				} else
					continue;
			}
			// if(RmeanGreatList.get(RmeanGreatList.size()-1)!=
			// meanGreatList.get(meanGreatList.size()-1))
			// {
			RmeanGreatList.add(meanGreatList.get(meanGreatList.size() - 1));
			// }
			// Combine two list RmeanLessList+ RmeanGreatList
			ArrayList<Integer> combineList = new ArrayList<>();
			combineList.addAll(RmeanLessList);
			combineList.addAll(RmeanGreatList);
			// add one more variable for extracting
			Collections.sort(combineList);
			combineList.add(ts.length+sIdx);
			// check combine_List
			/*
			 * ArrayList<Integer>combineList1=new ArrayList<>();
			 * combineList1.add(0); for(int i=1;i< combineList.size();i++){
			 * if((combineList
			 * .get(i-1)!=combineList.get(i))&&(combineList.get(i)-combineList
			 * .get(i-1))>2) combineList1.add(combineList.get(i)); }
			 */
			// display
			/*
			 * System.out.println("CombineList"); for(int i=0;i<
			 * combineList.size();i++){ System.out.print(combineList.get(i)+","); }
			 * System.out.println();
			 */
			// combineList1.remove(combineList1.size()-1);
			// combineList1.add(ts.length);
			/*
			 * System.out.println("CombineList1"); for(int i=0;i<
			 * combineList1.size();i++){ System.out.print(combineList1.get(i)+",");
			 * } System.out.println();
			 */
			return combineList;
			}	
		/*static class Sfeatures{
			double mean;
			double std;
			Sfeatures(double mean, double std){
				this.mean=mean;
				this.std=std;
			}
		}	*/	
		
	//PAA transform
	public static ArrayList<Double> PAA_Transform(double[] ts,ArrayList<Integer> SegArray, double nThreshold) {
		TSProcessor tsp = new TSProcessor();
		//ts = tsp.znorm(ts, nThreshold);
		ArrayList<Double> paa = new ArrayList<>();
		for (int i = 0; i < SegArray.size() - 2; i = i + 2) {
			double[] segment = Arrays.copyOfRange(ts, SegArray.get(i),SegArray.get(i + 2));
			paa.add(tsp.mean(segment));
		}
		return paa;
	}
		//Adaptive PAA distance
		  public static char[] SAX_Transform(double []ts, ArrayList<Integer>SegArray, double[]cuts,double nThreshold){
			  TSProcessor tsp=new TSProcessor();
			  /*//Display
			  for(int s=0;s<SegArray.size();s++){
				  System.out.print(SegArray.get(s)+" , ");
			  }
			  System.out.println();		  		 
			  ts=tsp.znorm(ts, nThreshold);		//normalization
			  ArrayList<Double>paa=new ArrayList<>();
			  for(int i=0;i< SegArray.size()-2;i=i+2){			 
					  double[]segment=Arrays.copyOfRange(ts, SegArray.get(i), SegArray.get(i+2));
					  paa.add(tsp.mean(segment));				  
				  }*/
			  ArrayList<Double>paa=PAA_Transform(ts,SegArray,nThreshold);
			  char[]sax_List=new char[paa.size()];
			  double []PAA=new double[paa.size()];
			  for(int i=0;i<paa.size();i++)
				  PAA[i]=paa.get(i);		  
			  sax_List=tsp.ts2String(PAA, cuts);
			  return sax_List;
		  }
	//adaptive SAX distance
		  public static double aSAX_distance(char[]ts, char[] qs, double[][]distanceMatrix, int n) throws SAXException{
			//check edge			
			  /*ArrayList<Integer>edge=new ArrayList<>();
				  for(int i=1;i<edgeList.size();i=i+2){
						  edge.add(edgeList.get(i)); 
					  }
				*/	  
			  if(ts.length== qs.length){
				  double dist=0.0D;
				  for(int i=1;i<ts.length;i++){				  
					  if(Character.isLetter(ts[i]) && Character.isLetter(qs[i])){
						  int numA=Character.getNumericValue(ts[i-1])-10;
						  int numB=Character.getNumericValue(qs[i-1])-10;
						  int maxIdx=distanceMatrix[0].length;
						  if(numA > (maxIdx-1) || numA < 0 ||numB > (maxIdx-1)||numB<0){
							  throw new SAXException("The Character index greater than "+ maxIdx + "or less than 0!");
						  }
						  double localDist=distanceMatrix[numA][numB];
						  //int tmp1=edge.get(i)-edge.get(i-1);
						  //dist=dist+ (tmp1*localDist*localDist);
						  dist=dist+ (localDist*localDist);
					  }
					  else{
						  throw new SAXException("Non-literal Character found!");
					  }
				  }
				  return Math.sqrt(dist);
			  }
			  else{
				  throw new SAXException("Data arrays lengths are not equal");
			  }
			  
		  }	

	// find Standard Deviation values for each paa segment
	public static double[] SD(double[] ts, ArrayList<Integer> SegArray) {
		TSProcessor tsp = new TSProcessor();
		ArrayList<Double> std = new ArrayList<>();
		for (int i = 0; i < SegArray.size() - 2; i = i + 2) {
			double[] segment = Arrays.copyOfRange(ts, SegArray.get(i),SegArray.get(i + 2));
			std.add(tsp.stDev(segment));
		}
		double[] stdArray = new double[std.size()];
		for (int i = 0; i < std.size(); i++) {
			stdArray[i] = std.get(i);
		}
		return stdArray;
	}
	
	// Compute symbol size using Minimum StandardDeviation in each TimeSeries  
				public static int MeanToResolutions(double maxVal) {
					double x_0 = 0; // minimum standard deviation
					double x_1 = 1.64; // maximum standard deviation
					int y_0 = 2; // minimum symbolic size
					int y_1 = 20; // maximum symbolic size
					int resolution = 0;
					double m=(y_1-y_0)/(double)(x_1-x_0);
					double temp_a1 = y_1-(Math.abs(maxVal)-x_0)*m;
						if ((int) Math.round(temp_a1) <= 0
								|| (int) Math.floor(temp_a1) == 1) 
							resolution = 2;
						else if ((int) Math.floor(temp_a1) >= 20)
							resolution = 20;
						else
							resolution = (int) Math.floor(temp_a1); // keep as absolute value of the alphabet size
					return resolution;
				}	
				  
	public static void main(String[] args) {
		if(args.length !=2 ){
			System.out.println("invalid argument: put 2-arguments for train and test file");
			System.exit(-1);
		}
		String train_filename=args[0];				
		String test_filename=args[1];
				
			long startTime=System.currentTimeMillis();
			List<sampleSeries>trainList=dataLoad(train_filename);
			List<sampleSeries>testList=dataLoad(test_filename);
			//SAXProcessor saxp=new SAXProcessor();
			TSProcessor tsp=new TSProcessor();
			Alphabet normalA=new NormalAlphabet();
			double nThreshold=0.0001;
			//int alpha=4;
			double bestsofar=Double.POSITIVE_INFINITY;
			int test_cLabel=-9999;
			int corrected=0;
			ArrayList<ArrayList<Integer>>sList=new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer>segList=new ArrayList<Integer>();
			double[][]dArray=new double[trainList.size()][trainList.get(0).Attributes.size()];
			double[][]tArray=new double[testList.size()][testList.get(0).Attributes.size()];
			double[]dummyArray=new double[trainList.get(0).Attributes.size()+1];
			for(int i=0;i<dummyArray.length;i++){
				dummyArray[i]=i*1D;			
			}
			//make Array
			for(int i=0;i<trainList.size();i++){
				for(int j=0;j< trainList.get(i).Attributes.size();j++){
					dArray[i][j]=trainList.get(i).Attributes.get(j);
				}
			}		
			
			for(int i=0;i<testList.size();i++){
				for(int j=0;j<testList.get(i).Attributes.size();j++){
					tArray[i][j]=testList.get(i).Attributes.get(j);
				}
			}		
			for(int i=0;i<dArray.length;i++){
				dArray[i]=tsp.znorm(dArray[i], nThreshold);
				segList=hmean_seg(dArray[i],0);
				sList.add(segList);
			}	
			
			//check Simple regressionLine
			ArrayList<ArrayList<Integer>> updatedSList = new ArrayList<ArrayList<Integer>>();
			for (int i = 0; i < sList.size(); i++) {
				/*
				 * for(int s=0;s<sList.get(0).size();s++){
				 * System.out.print(sList.get(0).get(s)+", "); }
				 */
				LinkedHashSet<Integer> updatedsegList = new LinkedHashSet<Integer>();

				double ethreshold = 0.05;
				for (int j = 0; j < sList.get(i).size() - 2; j = j + 2) {
					double[] tmp_y = Arrays.copyOfRange(dArray[i], sList.get(i)
							.get(j), sList.get(i).get(j + 2));
					double[] tmp_x = Arrays.copyOfRange(dummyArray, sList.get(i)
							.get(j), sList.get(i).get(j + 2));
					// display Array
					/*
					 * System.out.println("List"+ j); for(int
					 * a=0;a<tmp_y.length;a++){
					 * System.out.println(tmp_x[a]+" , "+tmp_y[a]); }
					 */
					SimpleRegression regression = new SimpleRegression();
					for (int s = 0; s < tmp_x.length; s++) {
						regression.addData(tmp_x[s], tmp_y[s]);
					}
					double mse = regression.getMeanSquareError();
					// System.out.println("Regression Mean Square Error:"+mse);
					if (mse > ethreshold) {
						updatedsegList.addAll(hmean_seg(tmp_y, sList.get(i).get(j)));

					} else {
						updatedsegList.add(sList.get(i).get(j));
						updatedsegList.add(sList.get(i).get(j + 1));
						updatedsegList.add(sList.get(i).get(j + 2));
					}

					regression.clear();
				}
				ArrayList<Integer> upsegList = new ArrayList<Integer>();
				upsegList.addAll(0, updatedsegList);
				updatedSList.add(upsegList);
			}
			Frequency q = new Frequency();	
			int tSeg=0;
			Set<Double>alphaSet=new HashSet<Double>();
			Map<Double,Integer>alphaMap=new HashMap<Double,Integer>();
			int maxCount=0;
			double maxVal=0D;
			//int minCount=Integer.MAX_VALUE;
			//double minVal=0D;
			int alphaSize=0;
			ArrayList<ArrayList<Double>> paa_List = new ArrayList<ArrayList<Double>>();
			for (int d = 0; d < dArray.length; d++)
			{
				ArrayList<Double>paa= PAA_Transform(dArray[d], updatedSList.get(d),nThreshold);
				paa_List.add(paa);
				tSeg += paa.size();
			}
			for (int i = 0; i < paa_List.size(); i++)
			{
				for (int j = 0; j < paa_List.get(i).size(); j++)
				{
					double tmp=Precision.round(paa_List.get(i).get(j), 2,BigDecimal.ROUND_HALF_UP);
					q.addValue(tmp);
					alphaSet.add(tmp);
				}
				for(Double s: alphaSet)
				{
					alphaMap.put(s, (int) q.getCount(s));
				}
				//display
				for(Map.Entry<Double, Integer>entry: alphaMap.entrySet()){
					/*if(minCount>entry.getValue()){
						minCount=entry.getValue();
						minVal=entry.getKey();
					}*/
					if(maxCount< entry.getValue()){
						maxCount=entry.getValue();
						maxVal=entry.getKey();
					}									
				
				}
			}
			/*double Gmean=0D;
			for(int i=0;i<dArray.length;i++){
				dArray[i]=tsp.znorm(dArray[i], nThreshold);
				for(int j=0;j<dArray[i].length;j++){
					Gmean+=dArray[i][j];
				}
			}*/
			alphaSize=MeanToResolutions(maxVal);		
			
			//classification
			for(int c=0;c< tArray.length;c++){
				bestsofar=Double.POSITIVE_INFINITY;
				for(int d=0;d< dArray.length;d++){
					//Symbolic Aggregate Approximation 				
					try {					
						double[]tSD_value=SD(dArray[d], updatedSList.get(d));
						double[]qSD_value=SD(tArray[c],updatedSList.get(d));
						char[] dSAX_List = SAX_Transform(dArray[d],updatedSList.get(d),normalA.getCuts(alphaSize),nThreshold);
						char[]tSAX_List=SAX_Transform(tArray[c],updatedSList.get(d),normalA.getCuts(alphaSize),nThreshold);
						
						double SDDist=0;
						for(int s=0;s<tSD_value.length;s++){
							//SDDist+=(double)(tSD_value.length/trainList.get(0).Attributes.size())*Math.pow((qSD_value[s]- tSD_value[s]), 2);
							SDDist+=Math.pow((qSD_value[s]- tSD_value[s]), 2);
						}
						double saxDist=aSAX_distance(dSAX_List,tSAX_List,normalA.getDistanceMatrix(alphaSize), dArray[d].length);
						double fsaxDist=Math.sqrt((double)(trainList.get(0).Attributes.size()/tSD_value.length))*(Math.sqrt(saxDist+SDDist));		
						if(fsaxDist< bestsofar){
							test_cLabel=trainList.get(d).cName;
							bestsofar=fsaxDist;
						}
					} catch (SAXException e) {
						e.printStackTrace();
					}				
				}
				if(test_cLabel==testList.get(c).cName)corrected=corrected+1;
			}
			System.out.println("Total execution time: "+(System.currentTimeMillis()-startTime) + " ms");	
			System.out.println("Alpha: " + alphaSize);
			System.out.println("total number of segments" + tSeg);
			System.out.println("Corrected Label: "+ corrected);
			System.out.println("The error rate is "+(double)(testList.size()-corrected)/(double)testList.size());
			System.out.println("***********************************************************");
		
	}
}	
