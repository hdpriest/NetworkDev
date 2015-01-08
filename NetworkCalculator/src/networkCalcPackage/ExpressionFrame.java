package networkCalcPackage;

import java.io.PrintWriter;

class ExpressionFrame {

    private int N;
    private int M;
    private float[][] DataFrame;
    private float[] k;
    private float[] means;
    private float[] gccSums;
    private String[] X_lab;
    private String[] Y_lab;
    private int X_iterator;

    public ExpressionFrame(int Dim1, int Dim2) {
        DataFrame = new float[Dim1][Dim2];
        N = Dim1;
        M = Dim2;
        k = new float[Dim1];
        means = new float[Dim1];
        gccSums = new float[Dim1];
        X_lab = new String[Dim1];
        Y_lab = new String[Dim2];
        X_iterator = -1;
    }

    public void resetIterator() {
        X_iterator = -1;
    }

    public int getNumColumns() {
        int I = DataFrame[0].length;
        return I;
    }

    public int getNumRows() {
        int I = DataFrame.length;
        return I;
    }

    public boolean hasNext() {
        if (DataFrame[X_iterator + 1] != null) {
            return true;
        } else {
            return false;
        }
    }

    public double[] getRowByIndexDbl(int I) {
        float[] Row = _getRow(I);
        double[] dRow = new double[M];
        for (int i = 0; i < M; i++) {
            dRow[i] = (double) Row[i];
        }
        return dRow;
    }

    private float[] _getRow(int I) {
        float[] Row = new float[M];
        for (int j = 0; j < M; j++) {
            Row[j] = _getValueByEntry(I, j);
        }
        return Row;
    }

    public float[] getRowByIndex(int I) {
        if (DataFrame[I] != null) {
            return _getRow(I);
        } else {
            System.err.println("Cannot get row " + I + " from matrix.\n\n");
            System.exit(0);
        }
        return null;
    }

    private float _getValueByEntry(int I, int J) {
        return DataFrame[I][J];
    }

    public float getValueByEntry(int I, int J) {
        return _getValueByEntry(I, J);
    }

    public void setValueByEntry(float Value, int I, int J) {
        _setValueByEntry(Value, I, J);
    }

    private void _setValueByEntry(float Value, int I, int J) {
        DataFrame[I][J] = Value;
    }

    public String[] getRowNames() {
        return X_lab;
    }

    public void setRowNames(String[] Rows) {
        X_lab = Rows;
    }

    public void setColumnNames(String[] Cols) {
        Y_lab = Cols;
    }

    public float[][] getDataFrame() {
        return DataFrame;
    }

    public float[] getNextRow() {
        float[] thisRow = new float[DataFrame[0].length];
        thisRow = DataFrame[X_iterator + 1];
        X_iterator++;
        return thisRow;
    }

    public void maskMatrix(float maskLevel) {
        int H = DataFrame.length;
        int W = DataFrame[0].length;
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                float v = _getValueByEntry(i, j);
                if (v < maskLevel) {
                    _setValueByEntry(0.0f, i, j);
                }
            }
        }
    }

    public boolean testValue(int i, int j) {
        boolean res = true;
        if (_getValueByEntry(i, j) == 0) {
            res = false;
        } else {
            res = true;
        }
        return res;
    }

    public void calculateKs() {
        int H = DataFrame.length;
        for (int i = 0; i < H; i++) {
            float thisK = 0f;
            for (int j = 0; j < DataFrame[i].length; j++) {
                if (i == j) {

                } else {
                    thisK += _getValueByEntry(i, j);
                }
            }
            k[i] = thisK;
        }
    }

    public float findK(int R, int j) {
        float K = k[R];
        return K;
    }
    private boolean _checkRow (float[] Data){
        // We preclude loading of all null/zero data rows
        // Also need to handle no-variance dataset (all 1's, all anything, this will throw many correlations for a loop)
        return true;
        // Added handling of these rows at similarity calculation
        /*
        float max = 0.0f;
        float min = Float.POSITIVE_INFINITY;
        //System.err.println("Checking Row");
        for(int i=0;i<Data.length;i++){
            if(Data[i] > max){
                max = Data[i];    
            }
            if(Data[i] < min){
                min = Data[i];
            }
            
        }
        if(max == min) return false; // occurrs on all-zero rows
        return true;
                */
    }
    public void addRow(float[] Row) {
        int I = X_iterator;
        if(_checkRow(Row)){
            System.arraycopy(Row, 0, DataFrame[I + 1], 0, Row.length);    
        }else{
            System.err.println("Row of null expression loaded. Handling for this occurrence is not implemented. Please add pseudocounts to null rows");
            System.err.println("Closing...");
            System.exit(0);
        }
        
        X_iterator++;
    }

    public void changeRow(int I, float[] Row) {
        System.arraycopy(Row, 0, DataFrame[I], 0, Row.length);
    }

    public void printMatrixToFile(String path, String Sep) {
        try {
            PrintWriter writer = new PrintWriter(path, "UTF-8");
            for (int y = 0; y < Y_lab.length; y++) {
                writer.print(Y_lab[y]);
                if (y != Y_lab.length - 1) {
                    writer.print(Sep);
                }
            }
            writer.print("\n");
            for (int i = 0; i < DataFrame.length; i++) {
                float[] Row = new float[DataFrame[i].length];
                Row = DataFrame[i];
                writer.print(X_lab[i]);
                writer.print(Sep);
                for (int j = 0; j < Row.length; j++) {
                    writer.print(Row[j]);
                    if (j != Row.length - 1) {
                        writer.print(Sep);
                    }
                }
                writer.print("\n");
            }
            writer.close();
        } catch (Exception e) {
            // 
        }

    }

    public void calculateMeans() {
        for (int i = 0; i < DataFrame.length; i++) {
            float[] Row = _getRow(i);
            float s = 0;
            for (int j = 0; j < Row.length; j++) {
                s += Row[j];
            }
            means[i] = s / Row.length;
        }

    }

    public float getMean(int i) {
        return means[i];
    }

    public float getGiniDenom(int i) {
        return gccSums[i];
    }

    public void calculateGiniSums() {
        for (int i = 0; i < DataFrame.length; i++) {
            float[] array1 = _getRow(i);
            float Denominator = 0.0f;
            int[] sortRanks1 = Operations.getIndicesInOrder(array1);
            for (int j = 0; j < array1.length; j++) {
                float v2 = ((2 * (j + 1)) - array1.length - 1) * array1[sortRanks1[j]];
                Denominator += v2;
            }
            gccSums[i] = Denominator;
        }
    }
}
