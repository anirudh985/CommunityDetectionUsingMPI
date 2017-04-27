package pc.vertex.code;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

public class RenumberNodes {
	
public static final int MAX_WORD_SIZE = 7;
	
	public static void main(String[] args) {
		VertexFollowing vf = new VertexFollowing();
		try {
			vf.readFile("C:\\Users\\Odin1\\Documents\\CommunityDetectionUsingMPI\\src\\out1.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void readFile(String filePath) throws IOException {
		FileReader fr = new FileReader(filePath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		
		Long counter = 0l;
		
		HashMap<Long, Long> numbersNew = new HashMap<Long, Long>();
		
		PrintWriter writer = new PrintWriter("C:\\Users\\Odin1\\Documents\\CommunityDetectionUsingMPI\\src\\out2.txt", "UTF-8");
		
		while((line = br.readLine()) != null) {
			String s[] = line.split("	");
			Long source = Long.valueOf(s[0]);
			Long destination = Long.valueOf(s[1]);
			if (numbersNew.containsKey(source)) {
	
				
			} else {
				numbersNew.put(source, ++counter);
			}
			
			writer.println(appendZero(numbersNew.get(source).toString(), 
					numbersNew.get(source).toString().length()) + "	" + 
					
					appendZero(destination.toString(), destination.toString().length()) +"	"+
					appendZeroForDouble(s[2].toString(), s[2].toString().length()));
			
		}
		
		
		FileReader fr2 = new FileReader(filePath);
		BufferedReader br2 = new BufferedReader(fr2);
		
		while((line = br2.readLine()) != null) {
			String s[] = line.split("	");
			Long source = Long.valueOf(s[0]);
			Long destination = Long.valueOf(s[1]);
			
			writer.println(appendZero(source.toString(), 
					source.toString().length()) + "	" + 
					appendZero(numbersNew.get(destination).toString(),
					numbersNew.get(destination).toString().length()) +"	"+
					appendZeroForDouble(s[2].toString(), s[2].toString().length()));
			
		}
		
		writer.close();
		br.close();
		br2.close();
	}
	
	private String appendZero(String s, int currLen) {
		for (int i = 0; i < MAX_WORD_SIZE - currLen; i++) {
			s = "0" + s;
		}
		return s;
	}

	private String appendZeroForDouble(String s, int currLen) {
		for (int i = 0; i < 2 * MAX_WORD_SIZE - currLen; i++) {
			s = "0" + s;
		}
		return s;
	}

}
