import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;

// this version Boyer Moore algorithm goes through the sequence, matching pattern characters from right to left side.
// part of preprocessing is generating a lookup map for the given pattern p
// when a mismatch occurs, the algorithm tries to skip some comparisons by:
// 1. Bad Character Rule - given a mismatched character s from sequence at pattern index i,
// the algorithm looks up the distance to the nearest s in pattern p using the lookup matrix
// which would be at i-th position of the row corresponding to character s
// Example: mismatched character is 'a' with index 5, distance to nearest 'a' to left side is
// map[5][0], since row [0] corresponds to 'a' and shifts the character comparison to right by
// that amount of characters
// 2. 'Good Character Rule' - given a mismatch after at least one character from pattern is matched,
// the algorithm looks up the distance to the first matched character (last character of pattern)
// then shifts the character comparison to right by that amount
//
// If both Bad Character Rule and Good Character Rule yield 0 as result, the algorithm
// shifts character comparisons by the length of pattern.

public class Main{

    static int[][] CharacterMap(String pattern) // generates a lookup map (implemented as integer matrix)
    // each value in map represents the distance to the nearest same character on left
    {
        int pLength = pattern.length();
        int[][] result = new int[pLength][5];

        for (int i = 0 ; i< pLength-1;i++) // last character does not change the map, thus is left out
        {
            switch (pattern.charAt(i))
            {
                case 'a': // [0], first row of matrix
                {
                    for (int q = i+1; q < pLength; q++)
                    {
                        result [q][0] = q - i;
                    }
                    break;
                }
                case 'c': // [1], second row of matrix
                {
                    for (int q = i+1; q < pLength; q++)
                    {
                        result [q][1] = q - i;
                    }
                    break;
                }
                case 't': // [2], third row of matrix
                {
                    for (int q = i+1; q < pLength; q++)
                    {
                        result [q][2] = q - i;
                    }
                    break;
                }
                case 'g': // [3], fourth row of matrix
                {
                    for (int q = i+1; q < pLength; q++)
                    {
                        result [q][3] = q - i;
                    }
                    break;
                }
                case 'n': // [4], fifth row of matrix
                {
                    for (int q = i+1; q < pLength; q++)
                    {
                        result [q][4] = q - i;
                    }
                    break;
                }
                default:
                    System.out.print("Character not detected\n");
                    break;
            }
        }

        return result;
    }

    static int CharacterIndex(char z, int[][] map, int patternIndex) // used to look up Char in preprocessed map
    {
        // looks up how much distance to the same character on the right using lookup map
        switch (z)
        {
            case 'a': // in first row
            {
                return map[patternIndex][0];
            }
            case 'c': // in second row
            {
                return map[patternIndex][1];
            }
            case 't': // in third row
            {
                return map[patternIndex][2];
            }
            case 'g': // in fourth row
            {
                return map[patternIndex][3];
            }
            case 'n': // fifth row
            {
                return map[patternIndex][4];
            }
            default:
                break;
        }
        return 0;
    }



    static private int[] boyerMoore(String sequence, String pattern)
    {

        int patternPositions[] = new int[10];
        int map[][], sLength = sequence.length(),pLength=pattern.length();
        int i,j,bcr,gcr,shiftL,matchCounter=0;
        boolean ms; //ms for mismatch
        int shiftCounter = 0;
        String temp="", temp2 = "";
        map = CharacterMap(pattern); // generate character map for current pattern
        for (i = 0 ; i < sLength-pLength-1; i++)
        {
            ms = false;
            shiftL = -1;
            temp2 = "";
            for (j = pLength-1; j>=0; j--)
            {
                if (pattern.charAt(j)!=sequence.charAt(i+j))
                {
                    ms = true;
                    bcr = CharacterIndex(sequence.charAt(i+j),map,j); // bad character rule
                    gcr = CharacterIndex(pattern.charAt(pLength-1),map,pLength-1); // good character rule
                    // instead of good suffix rule, we implemented a 'Good Character Rule' that
                    // only takes the last character of pattern as suffix and sets gcr equal
                    // to the amount of characters that can be skipped according to lookup map

                    shiftL = bcr;
                    if (gcr>bcr && j < pLength-1){ // set shift equal to the greater number between bcr and gcr
                            shiftL = gcr;
                    }
                    if (gcr == 0 && bcr ==0) // if character is not in pattern, skip next pLength positions, (-1 because i iterates)
                    {
                        i+= pLength-1;
                        break;
                    }
                    if (shiftL>0) { // slide matching index 'i' to the right by 'Sh
                        if (i+shiftL > sLength - 1) // end the algorithm if i+shiftL is past the length of the sequence
                        {return patternPositions;}
                        shiftCounter += shiftL-1;
                        i += shiftL - 1;
                        break;
                    }// -1 because i is going to increment by 1 anyway*/
                }
            }
            if (!ms) // if all characters match, then we found the pattern in sequence
            {
                if (matchCounter > 9) { // increase the total number of matches found
                    matchCounter++;
                } else {
                    patternPositions[matchCounter++] = i+1; // store first 10 values in patternPositions
                }
            }
        }
        //System.out.print(shiftCounter); // amount of comparisons skipped
        System.out.print(" " + matchCounter + '\n');  // total number of matches
        return patternPositions;
    }

    public static void main(String[] args)
    {
        int[] results = new int[10];
        try
        {
            String temp = "";
            // import sequence
            String sequence = "";
            BufferedReader sequenceReader = new BufferedReader(new FileReader(args[1]));
            StringBuilder importedSequence = new StringBuilder();
            while ((temp=sequenceReader.readLine())!=null)
            {
                if (temp.charAt(0)!='>')
                {
                    importedSequence.append(temp);
                }
                temp = "";
            }
            sequence = importedSequence.toString();

            // import patterns
            BufferedReader patternReader = new BufferedReader(new FileReader(args[0]));
            StringBuilder importedPatterns = new StringBuilder();
            long timeSpent = System.currentTimeMillis();
            while (patternReader.readLine()!=null)
            {
                temp=patternReader.readLine();
                if (temp.charAt(0)!='>')
                {

                    System.out.print(temp + ": ");
                    results = boyerMoore(sequence,temp); // call boyerMoore method on current pattern and the sequence
                }

                System.out.print("[");
                for (int p = 0 ; p < 10; p++)
                {
                    System.out.print(results[p]); // print stored indices
                    if (p<9) {System.out.print(", ");}
                }
                System.out.print("]\n");
                temp = "";
            }
            timeSpent = System.currentTimeMillis() - timeSpent; // total time for matching all patterns (in milliseconds)
            System.out.print(timeSpent + " milliseconds");
            patternReader.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        catch (Exception E) {
            System.err.println(E);
            E.printStackTrace();
        }
    }


}
