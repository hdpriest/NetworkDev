package networkCalcPackage;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.ujmp.core.enums.FileFormat;
import org.ujmp.core.exceptions.MatrixException;
import org.ujmp.core.floatmatrix.FloatMatrix2D;
import org.ujmp.core.interfaces.GUIObject;
import org.ujmp.gui.io.ExportPNG;




class GCNMatrix {
	private int N;
	private float[][] DataFrame;
	private float[] k;
	private float[] means;
	private float[] gccSums;
	private String[] X_lab;
	private String[] Y_lab;
	private int X_iterator;
	
	
         public void generateHeatmap () {
           double[][] newDF = new double[N][N];
           FloatMatrix2D nat = FloatMatrix2D.factory.dense(N, N);       
           for(int i=0;i<N;i++){
               float[] Row = getRowByIndex(i);
               for(int j=0;j<Row.length;j++){
                   nat.setFloat(Row[j],i,j);
               }
           }
           GUIObject gui = nat.getGUIObject();
          JFrame frame = gui.getFrame();
          
          BufferedImage img = new BufferedImage(gui.getFrame().getWidth(),gui.getFrame().getHeight(),BufferedImage.TYPE_INT_RGB);
          Graphics2D g2d = img.createGraphics();
          File outputfile = new File("image.png");
            try {
                ImageIO.write(img, "png", outputfile);
            } catch (IOException ex) {
                Logger.getLogger(GCNMatrix.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
	public GCNMatrix (int Dim1, int Dim2) {
		DataFrame = new float[Dim1][Dim2];
		N = Dim1;
		k = new float[Dim1];
		means = new float[Dim1];
		gccSums= new float[Dim1];
		X_lab = new String[Dim1];
		Y_lab = new String[Dim2];
		X_iterator = -1;
	}
	
	public double[] generateHistogram (String pathOut, String Title,String Xlab, String Ylab) {
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		double[] Histogram=new double[201];
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.HALF_UP);
		for(int i=0;i<H;i++){
			for(int j=i;j<W;j++){
				//System.out.println("Val: "+DataFrame[i][j]+"\n");
				if(DataFrame[i][j] != 0.0){
					Double v= ((Double.valueOf(df.format(DataFrame[i][j])))+1)*100;
					int value = v.intValue();
					Histogram[value]++;
				}
			}
		}
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Values");
		for(int x=0;x<201;x++){
			double a = ((double)x/100)-1.0;
			Double A = Double.valueOf(df.format(a));
			//System.out.println("Index: " + x + " Value: "+ A +"\tObs: "+Histogram[x]);
			series.add((double) A, (double) Histogram[x],true);			
		}
		dataset.addSeries(series);
		JFreeChart chart = ChartFactory.createXYLineChart(Title,Xlab,Ylab,dataset, PlotOrientation.VERTICAL, 
				 false, true, false);
		final Plot plot = chart.getPlot();
		plot.setBackgroundPaint(Color.white);
		plot.setOutlinePaint(Color.black);
		chart.setBackgroundPaint(Color.white);
		chart.setAntiAlias(true);
		try {
			ChartUtilities.saveChartAsJPEG(new File(pathOut), chart, 500, 300);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		return Histogram;
	}
	
	public void resetIterator () {
		X_iterator=-1;
	}
	
	public int getNumColumns () {
		int I = DataFrame[0].length;
		return I;
	}
	
	public int getNumRows () {
		int I = DataFrame.length;
		return I;
	}
	
	public boolean hasNext () {
		if(DataFrame[X_iterator+1] != null){
			return true;
		}else{
			return false;			
		}
	}
	public double[] getRowByIndexDbl (int I){
		float[] Row = _getRow(I);
		double[] dRow = new double[N];
		for(int i=0;i<Row.length;i++){
			dRow[i] = (double) Row[i];
		}
		return dRow;
	}
        
	private float[] _getRow(int I){
		float[] Row = new float[N];
		for(int j=0;j<N;j++){
			Row[j]=_getValueByEntry(I,j);
		}
		return Row;
	}
        
        private float[] _getRowAsDistance(int I){
		float[] Row = new float[N];
		for(int j=0;j<N;j++){
			Row[j]=1.0f - (_getValueByEntry(I,j));
		}
		return Row;
	}
        
	public float[] getRowByIndexAsDistance (int I){
		if(DataFrame[I] != null){
			return _getRowAsDistance(I);
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return null;
	}
        
	public float[] getRowByIndex (int I){
		if(DataFrame[I] != null){
			return _getRow(I);
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return null;
	}
	
	private float _getValueByEntry (int I,int J){
		if(I<=J){
			return DataFrame[I][J];
		}else{
			return DataFrame[J][I];
		}
	}
        
        private float _getValueByEntryAsDistance (int I,int J){
		if(I<=J){
			return (1.0f-DataFrame[I][J]);
		}else{
			return (1.0f-DataFrame[J][I]);
		}
	}
	
	public float getValueByEntry (int I,int J){
		return _getValueByEntry(I,J);
	}
        
        public float getValueByEntryAsDistance (int I,int J){
		return _getValueByEntryAsDistance(I,J);
	}
	
	public void setValueByEntry (float Value,int I, int J){
		_setValueByEntry(Value,I,J);
	}
	
	private void _setValueByEntry (float Value,int I, int J){
		if(I<=J){
			DataFrame[I][J]=Value;
		}else{
			DataFrame[J][I]=Value;
		}
	}
	
	public String[] getRowNames (){
            return X_lab;
        }
	public String getRowName (int I){
        return X_lab[I];
    }
	public void setRowNames (String[] Rows) {
            X_lab = Rows;
	}
	
	public void setColumnNames (String[] Cols){
		Y_lab = Cols;
	}
	public float[][] getDataFrame () {
            return DataFrame;
        }
        
	public float[] getNextRow () {
		float[] thisRow = new float[DataFrame[0].length];
		thisRow=DataFrame[X_iterator+1];
		X_iterator++;
		return thisRow;
	}
	
	public void maskMatrix (float maskLevel) {
		int H = N;
		for(int i=0;i<H;i++){
			for(int j=i;j<N;j++){
				float v = _getValueByEntry(i,j);
				if(Math.abs(v)<maskLevel){
                                        _setValueByEntry(0.0f,i,j);
				}
			}
		}
	}
        
       
        
	public TreeMap<Float,Integer> generateDistribution () {
	int H = N;
	DecimalFormat df = new DecimalFormat("#.##");
	df.setRoundingMode(RoundingMode.HALF_UP);
	TreeMap<Float,Integer> HMHistogram = new TreeMap<Float,Integer>();
	for(int i=0;i<H;i++){
		for(int j=i;j<H;j++){
			//System.out.println("Val: "+DataFrame[i][j]+"\n");
			if(_getValueByEntry(i,j) != 0){
				try{
					Float V = (Float.valueOf(df.format(_getValueByEntry(i,j))));
					if(HMHistogram.containsKey(V)){
						Integer I = HMHistogram.get(V);
                                               // System.out.println("Putting " + V+ " and " + I);
						HMHistogram.put(V,I+1);
					}else{
						HMHistogram.put(V,1);
                                               // System.out.println("Putting " + V+ " and 1");
					}
				}catch(NumberFormatException ex){
					System.out.println("Obtain " + _getValueByEntry(i,j) +" from matrix.");
					System.exit(1);
				}
			}
		}
	}
	return HMHistogram;
}
  
	public boolean testValue (int i,int j){
		boolean res=true;
		if(_getValueByEntry(i,j) == 0f){
			res=false;
		}else{
			res=true;
		}
		return res;
	}
	
	public void calculateKs (){
		int H = DataFrame.length;
		for(int i=0;i<H;i++){
			float thisK=0f;
			for(int j=0;j<DataFrame[i].length;j++){
				if(i==j){
					
				}else{
					thisK+=_getValueByEntry(i,j);
				}
			}
			k[i]=thisK;
		}	
	}
	
	public float findK (int R,int j){
		/*for(int i=0;i<DataFrame[R].length;i++){
			if(i==j){
				/// do not add DataFrame[R][j] to the k of R
			}else{
				K+=DataFrame[R][i];
			}
		}
		return K;*/
		float K=k[R];
		return K;
	}
	
	public void addRow (float[] Row){
		int I=X_iterator;
                System.arraycopy(Row, 0, DataFrame[I+1], 0, Row.length);
		X_iterator++;
	}
	
	public void changeRow (int I,float[] Row){
                System.arraycopy(Row, 0, DataFrame[I], 0, Row.length);
	}
	
	public void printMatrixToFile (String path,String Sep){
		try {
			PrintWriter writer = new PrintWriter(path,"UTF-8");
			writer.print(Sep);
			for(int y=0;y<Y_lab.length;y++){
				writer.print(Y_lab[y]);
				if(y != Y_lab.length-1){
					writer.print(Sep);
				}
			}
			writer.print("\n");
			for(int i=0;i<DataFrame.length;i++){
				float[] Row = new float[DataFrame[i].length];
				Row = DataFrame[i];
				writer.print(X_lab[i]);
				writer.print(Sep);
				for(int j=0;j<Row.length;j++){
					writer.print(Row[j]);
					if(j != Row.length-1){
						writer.print(Sep);
					}
				}
				writer.print("\n");
			}
			writer.close();
		} catch (Exception e){
			// 
		}
		
	}

	public void calculateMeans() {
		for(int i=0;i<DataFrame.length;i++){
			float[] Row = _getRow(i);
			float s=0;
			for(int j=0;j<Row.length;j++){
				s+=Row[j];
			}
			means[i]=s/Row.length;
		}
		
	}

	public float getMean(int i){
		return means[i];
	}
	
	public float getGiniDenom(int i) {
		return gccSums[i];
	}
	
	public void calculateGiniSums() {
		for(int i=0;i<DataFrame.length;i++){
			float[] array1 = _getRow(i);
			float Denominator=0.0f;
			int[] sortRanks1 = Operations.getIndicesInOrder(array1);
			for(int j=0;j<array1.length;j++){
				float v2=((2*(j+1))-array1.length-1) * array1[sortRanks1[j]];
				Denominator+=v2;
			}
			gccSums[i]=Denominator;
		}
	}
}
