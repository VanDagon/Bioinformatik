import java.util.*;
import java.io.File;
import java.io.*;

public class Clustal {

    public static int[][] generateInitialMatrix(int nSequences, String[] importedSequences)
    {
        int counter = 1;
        String[] results;
        globalAlignment GA = new globalAlignment();
        int[][] matrix = new int[nSequences][nSequences];
        for (int i = 0 ; i < nSequences; i ++)
        {
            for (int j = 0; j < nSequences ; j ++)
            {
                results = GA.findAlignment(importedSequences[i],importedSequences[j]);
                matrix[i][j] = GA.bestscore * (i==j ? 0:1);
                GA.getMRID(results);
                GA.initializeGA();
            }
            counter++;
        }
        return matrix;
    }

    public static int[] findMax(int[][] matrix)
    {
        int max = matrix[0][0];
        int[] results = new int[2];
        for (int i=0 ; i<matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++){
                if (matrix[i][j] > max) {max = matrix[i][j]; results[0] = i; results[1] = j;}
            }
        }
        return results;
    }

    public static int[][] hierarchicalClustering( int nSequences, String[] importedSequences, String fileName)
    {
        boolean fileStart = false;
        String[] names = globalAlignment.readNames(fileName,8);
        int[][] m = generateInitialMatrix(nSequences,importedSequences);
        printMatrix(m,names);
        ArrayList<String> graph = new ArrayList<String>();
        for (int i = 0 ; i < names.length; i ++)
        {
            graph.add(names[i]);
        }

        int ml = m[0].length;
        int minMax,maxMax;
        int[] maxCoor = findMax(m);
        for (int counter = 1; counter < ml -1 ; counter ++ )
        {
            int[][] temp = new int[ml-counter][ml-counter];
            minMax = Integer.min(maxCoor[0],maxCoor[1]);
            maxMax = Integer.max(maxCoor[0],maxCoor[1]);
            for (int i = 0 ; i < ml-counter; i ++) {
                for (int j = 0; j < ml-counter; j++) {
                    if (i==j) {continue;}

                    if (i == minMax && j!=maxMax) {
                        temp[i][j] = (m[i][j] + m[maxMax][j])/2;
                    } else if (j == minMax && i!=maxMax) {
                        temp[i][j] = (m[i][j] + m[i][maxMax])/2;
                    } else if (i >= maxMax && j>= maxMax) {
                        temp[i][j] = m[i+1][j+1];
                    } else if (j >= maxMax){
                        temp[i][j] = m[i][j+1];
                    } else if (i >= maxMax){
                        temp[i][j] = m[i+1][j];
                    } else {
                        temp[i][j] = m[i][j];
                    }
                }

            }
            m = new int[ml-counter][ml-counter];
            for (int q = 0 ; q < m.length; q++)
            {
                for (int z =0 ; z < m[0].length; z++)
                {
                    m[q][z] =temp[q][z];
                }
            }
            maxCoor = findMax(temp);
            System.out.println();
            try {
                String[] merged = new String[3];
                merged[0] = graph.get(minMax)+graph.get(maxMax);
                merged[1] = graph.get(minMax);
                merged[2] = graph.get(maxMax);
                File file = new File("output.txt");
                BufferedWriter bw = new BufferedWriter(new FileWriter(file, true));
                if (!fileStart){fileStart=true;bw.append("Digraph G {");bw.newLine();}
                bw.append(merged[0] + " - > " + merged[1]);
                bw.newLine();
                bw.append(merged[0] + " - > " + merged[2]);
                bw.newLine();
                bw.close();
            }
         catch (IOException e) {
        e.printStackTrace();
        } catch (Exception E) {
        System.err.println(E);
        E.printStackTrace();
        }
            graph.set(minMax,(graph.get(minMax)+graph.get(maxMax)));
            graph.remove(graph.get(maxMax));
        }

        try {
            File file = new File("output.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(file, true));
            if (!fileStart){fileStart=true;bw.append("Digraph G {");bw.newLine();}
            bw.append(graph.get(0)+graph.get(1) + " - > " + graph.get(0));
            bw.newLine();
            bw.append(graph.get(0)+graph.get(1) + " - > " + graph.get(1));
            bw.newLine();
            bw.append("}");
            bw.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        } catch (Exception E) {
            System.err.println(E);
            E.printStackTrace();
        }
        return m;
    }

    public static void printMatrix(int[][] matrix, String[] names)
    {
        int xl = names.length;
        int maxL = names[0].length();
        for (int i = 0 ; i < xl ; i++)
        {
            if (names[i].length()>maxL) {maxL = names[i].length();}
        }
        for (int i = 0 ; i < maxL; i++) {
            System.out.print(' ');
        }
        for (int i = 0 ; i < xl ; i++)
        {
            do {names[i]+=' ';}
            while (names[i].length()<maxL+1);
            System.out.print(names[i]);
        }
        System.out.println();
        int temp, counter;
        for (int i=0; i<xl;i++)
        {
            System.out.print(names[i]);
            for (int j = 0 ; j < xl; j++)
            {
                counter = 1;
                temp = matrix[i][j];
                System.out.print(matrix[i][j]);
                while (temp>10){temp/=10;counter++;}
                for (int q = 0 ; q < maxL-counter+1; q++) {
                    System.out.print(" ");
                }
            }
            System.out.println();

        }
    }


    public static void main(String[] args) {
        int nSequences = 8;
        globalAlignment GA = new globalAlignment();
        String[] importedSequences = GA.readFasta(args[0],nSequences); // "src/pair.fasta"
        GA.readScoringMatrix(args[1]); // "src/matrix.txt"
        String[] names = GA.readNames(args[0],nSequences); // "src/pair.fasta"

        hierarchicalClustering(nSequences,importedSequences,args[0]);
    }
}
