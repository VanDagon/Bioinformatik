import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;

public class globalAlignment {
    private static int delCost = -4; // deletion cost, can be adjusted with args[2]
    private static int extendCost = -1; // not implemented
    private static String xchars="",ychars=""; //  xchars is list of chars in scoring  matrix horizontally, ychars vertically
    public static int[][] scoreMatrix; // = {{5, -3, -4, -4}, {-3, 5, -4, -4}, {-4, -4, 5, -2}, {-4, -4, -2, 5}};
    private static int blal = 0 ; //  short for Best Local Alignment Length
    private static String indexString = ""; // string that stores indices of best local alignments (if there are multiple with same score)
    public static int bestscore = 1; // s
    private static int nomat, norep, nodel, noins; // matches, replacements, deletions, insertion

    public void initializeGA() //initializes variables of global alignment
    {
        blal = 0;
        bestscore = 1;
        nomat = 0;
        norep = 0;
        nodel = 0;
        noins = 0;
    }

    public static String[] readNames(String fileName,int sequenceCount) {
        String[] names = new String[sequenceCount];
        int tempi,tempj;
        for (int i = 0; i < sequenceCount; i++)
        {
            names[i] = "";
        }
        try {
            String temp = "";
            BufferedReader sequenceReader = new BufferedReader(new FileReader(fileName));
            int counter = -1;
            while ((temp = sequenceReader.readLine()) != null) {
                if (temp.charAt(0) == '>') {
                    counter++;
                    tempi = temp.indexOf('[')+1;
                    tempj = temp.indexOf(' ',tempi);
                    names[counter] = temp.substring(tempi,tempj);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception E) {
            System.err.println(E);
            E.printStackTrace();
        }
        return names;
    }

    public static String[] readFasta(String fileName,int sequenceCount) {
        String[] importedSequence = new String[sequenceCount];
        for (int i = 0; i < sequenceCount; i++)
        {
            importedSequence[i] = "";
        }
        try {
            String temp = "";
            BufferedReader sequenceReader = new BufferedReader(new FileReader(fileName));
            int counter = -1;
            while ((temp = sequenceReader.readLine()) != null) {
                if (temp.charAt(0) != '>') {
                    importedSequence[counter] += temp;
                }
                else
                {
                    counter++;
                    temp = "";
                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception E) {
            System.err.println(E);
            E.printStackTrace();
        }
        return importedSequence;
    }

    public static int[][] readScoringMatrix(String fileName)
    {
        /* ignores lines that begin with #
           skips the lines that are empty
           expects list of characters at the first non-empty line
           expects n amount of following lines (where n is the length of character list)
           expects each following line to start with a single char
           n amount of values follow, corresponding to the tuple of chars
         */

        int[][] defaultResult = new int[1][1]; // placeholder matrix, just to satisfy a return for method
        defaultResult[0][0] = -99; // this should not be returned

        boolean charsRead = false; // true if the first line of matrix (chars corresponding to matrix columns) has been read
        int noc=0; // number of characters
        String temp=""; // temp stores readLine

        try {
            BufferedReader sequenceReader = new BufferedReader(new FileReader(fileName));

            while (!charsRead) // terminate if the line containing only chars has already been read
            {
                int xi = 0;
                temp = sequenceReader.readLine();
                if (temp.length()==0){continue;}
                if (temp.charAt(0)!= '#')
                {
                    while (xi<temp.length())
                    {
                        if (temp.charAt(xi)!= ' ')
                        {
                            if (!charsRead) {charsRead = true;}
                            xchars+=temp.charAt(xi); // append the found char to xchars
                        }
                        xi++;
                    }
                }
            }

            noc = xchars.length(); // number of chars is the length of xchars
            /* initialize matrix with number of columns = number of rows = number of chars in xchars */
            int result[][] = new int[noc][noc]; // who's there?
            boolean charRead; // true if the char corresponding to the matrix row has been read for each line
            temp = ""; // stores readLine
            int xi, counter = 0; // xi used to iterate through temp, counter is for rows of matrix
            int valCounter, nextSpace; // for column of matrix, nextSpace for index of next space

            while ((temp = sequenceReader.readLine()) != null && temp.length()!=0 )// terminate if there are no more lines & this line is empty
            {
                charRead = false;
                xi = 0;
                valCounter = 0;

                while (xi<temp.length() || valCounter!=noc) // terminate when either reached the end of line, or found enough values
                {
                    if (temp.charAt(xi)!=' ') // non-space character has been found
                    {
                        if (!charRead) { ychars+= temp.charAt(xi);charRead = true;}
                        else {
                            nextSpace = temp.indexOf(' ',xi);
                            if (nextSpace == -1) {nextSpace = temp.length();} // if no space is found, then we are at the end of line
                            /* takes substring from charAt(xi) until the next space, parses it to int and writes in matrix */
                            result[counter][valCounter] = Integer.parseInt(temp.substring(xi,nextSpace));
                            xi = temp.indexOf(' ',xi); // skip some comparisons
                            if (xi<0){break;}
                            valCounter++;
                        }
                    }
                    xi++;
                }
                counter++;
            }
            scoreMatrix = result;
            return result;
        }
        catch (IOException e) {
            e.printStackTrace();
        }

        return defaultResult;
    }

    public static int lookUpCost(char a, char b) {

        int al = xchars.length()-1, bl = ychars.length()-1;
        int ai = -1, bi = -1; // ai is column index of character a, bi is row index of character b
        for (int i = 0 ; i < xchars.length(); i ++)
        {
            if (xchars.charAt(i)==a){ai = i;}
        }
        for (int i = 0 ; i < ychars.length(); i ++)
        {
            if (ychars.charAt(i)==b){bi = i;}
        }
        // if neither of characters were found, give back the values of * and *
        if (ai < 0 && bi < 0) {
            return scoreMatrix[al][bl];
        }
        if (ai < 0) { //
            return scoreMatrix[bi][al];
        }
        if ( bi < 0) {
            return scoreMatrix[bl][ai];
        }
        return scoreMatrix[bi][ai];
    }

    static int getMax(int i, int j, int[][] aMatrix, String a, String b) // aMatrix --> alignment matrix (not fully generated)
    {
        int left, up, leftup, max; // value derived from left, up, or leftup (diagonal)
        if (i > a.length() || j > b.length()) {
            System.out.println("Wrong parameters");
            return -1;
        }
        left = aMatrix[i][j - 1] + delCost; // deletion/insertion from left
        up = aMatrix[i - 1][j] + delCost; // deletion/insertion from up
        leftup = aMatrix[i - 1][j - 1] + lookUpCost(a.charAt(i), b.charAt(j)); // match/mismatch from leftup (diagonal)



        max = Integer.max(left, Integer.max(up,leftup)); // take the max from the three values
        if (max > bestscore) { // if this is new best local alignment score
            bestscore = max; // overwrite bestscore
            /* empty indexString and write current indices in it */
            indexString = Integer.toString(i).concat(",").concat(Integer.toString(j)).concat(";");
            return max;
        }
        if (max == bestscore) { // if this score is the same as previous best, append these indices to indexString
            indexString.concat(Integer.toString(i)).concat(",").concat(Integer.toString(j)).concat(";");
        }
        return max;
    }


    static int[][] alignmentMatrix(String a, String b) // generates the alignment matrix
    {
        /* definitions / initializations */
        int xl = a.length() , yl = b.length() ; // get string lengths (+1)
        int[][] matrix = new int[xl][yl]; // initialize matrix
        matrix[0][0] = 0;
        /* generating first row and column of matrix */
        for (int i = 1; i < yl; i++) // generate insertions for first row
        {
            matrix[0][i] = delCost * i; // note that delCost stands for insertion AND deletion cost
        }

        for (int i = 1; i < xl; i++) // generate insertions for first column
        {
            matrix[i][0] = delCost * i; // note that delCost stands for insertion AND deletion cost
        }

        /* matrix generation */
        for (int i = 1; i < xl; i++) // begin generation horizontally
        {

            for (int j = 1; j < yl; j++) {
                /* look up the cost corresponding to the character tuple */
                matrix[i][j] = getMax(i, j, matrix, a, b);
            }
        }

        return matrix;
    }

    public static void printAlignment(String[] alignedSequences) {
        for (int i = 0; i < alignedSequences.length-1; i++) {
            for (int k = 0; k+50 < alignedSequences[i].length(); k+=50) {
                System.out.println(alignedSequences[i].substring(k, k+50));
                //print the signs for match/insert/replace
                for (int j = k; j < k+50; j++) {
                    if (alignedSequences[i].charAt(j) == alignedSequences[i + 1].charAt(j)) {
                        System.out.print("|");
                        nomat++;
                    } else if ((alignedSequences[i].charAt(j) == '_') || (alignedSequences[i + 1].charAt(j) == '_')) {
                        System.out.print(" ");
                        noins++;
                    } else {
                        System.out.print(".");
                        norep++;
                    }
                }
                System.out.println();
                System.out.println(alignedSequences[i+1].substring(k, k+50));
            }
            int l = alignedSequences[i].length()-alignedSequences[i].length()%50-1;
            System.out.println(alignedSequences[i].substring(l));
            for(int q=l; q < alignedSequences[i].length(); q++){
                if (alignedSequences[i].charAt(q) == alignedSequences[i + 1].charAt(q)) {
                    System.out.print("|");
                    nomat++;
                } else if (alignedSequences[i].charAt(q) == '_') {
                    System.out.print(" ");
                    noins++;
                } else if(alignedSequences[i+1].charAt(q) == '_'){ System.out.print(" ");
                    nodel++;}
                else {
                    System.out.print(".");
                    norep++;
                }
            }
            System.out.println();
            System.out.println(alignedSequences[i+1].substring(l));
        }
    }

    public static void getMRID(String[] alignedSequences) {
        for (int i = 0; i < alignedSequences.length-1; i++) {
            for (int j = 0; j < alignedSequences[i].length() - 1; j++)
            {
                if (alignedSequences[i].charAt(j) == alignedSequences[i + 1].charAt(j)) {
                    nomat++;
                } else if (alignedSequences[i].charAt(j) == '_') {
                    noins++;
                } else if (alignedSequences[i + 1].charAt(j) == '_'){
                    nodel++;
                } else {
                    norep++;
                }
            }
        }
    }

    public static String append(char c, String s) // append String s to char c ( result = c + s)
    {
        String temp;
        temp = c + s;
        return temp;
    }

    public static void printStats() {
        System.out.println("Length of best alignment: " +blal);
        System.out.println("Score of best alignment: " +bestscore);
        System.out.println("Number of matches: " +nomat);
        System.out.println("Number of Replacements: " +norep);
        System.out.println("Number of insertions: " +noins);
        System.out.println("Number of deletsions: " +nodel);
    }

    public static String[] findAlignment(String a, String b) { // traceback fuction for finding alignment
        String[] results = new String[2];
        results[0] = "";
        results[1] = "";
        int[][] aMatrix = alignmentMatrix(a, b);
        String tempi, tempj;
        int iindex=a.length()-1, jindex=b.length()-1;
        int l=1;

            /* aMatrix[iindex][jindex] is one of the maximum local alignment scores */
            while (iindex != 1  && jindex != 1) // start moving backwards from iindex,jindex
            {
                /* store the maximum value of left, leftup, up in max */
                int max = Integer.max(aMatrix[iindex - 1][jindex], Integer.max(aMatrix[iindex - 1][jindex - 1], aMatrix[iindex - 1][jindex - 1]));
                if (aMatrix[iindex - 1][jindex - 1] == max) // if leftup is max value
                {
                    results[0] = append(a.charAt(iindex - 1), results[0]);
                    results[1] = append(b.charAt(jindex - 1), results[1]);
                    iindex--;
                    jindex--;
                } else if (aMatrix[iindex - 1][jindex] == max) // if left is max value
                {
                    results[0] = append(a.charAt(iindex - 1), results[0]);
                    results[1] = append('_', results[1]);
                    iindex--;
                } else if (aMatrix[iindex][jindex - 1] == max) // if up is max value
                {
                    results[0] = append('_', results[0]);
                    results[1] = append(b.charAt(jindex - 1), results[1]);
                    jindex--;
                }
            }

        return results;
    }


/*    public static void main(String[] args) {
        scoreMatrix = readScoringMatrix("src/new.txt");
        String[] alignedSequences = readFasta("src/pair.fasta",2);
        String[] r = findAlignment(alignedSequences[0],alignedSequences[1]);
        printAlignment(r);


    }*/
}
