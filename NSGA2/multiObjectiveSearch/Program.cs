
using System;
using System.Collections.Generic;
using System.IO;

namespace multiObjectiveSearch
{
	class Program
	{
		public static void Main(string[] args)
		{
			
			NSGA2 n = new NSGA2("settings.txt");
			List<chromosome> ansn = n.SearchDesignSpace();
            
			StreamWriter sw = null;
			
			try
			{
				sw = new StreamWriter("NSGA2_output.txt");
			}
			catch
			{
				if(sw != null) 
					sw.Close();
			}
			if(sw != null)
			{
				for(int i = 0; i < ansn.Count; i++)
				{
					//sw.WriteLine("ans " + i.ToString() + " : " + ans[i].PrintCMDString());
					sw.WriteLine(ansn[i].PrintRawString("\t"));
				}
			}
			sw.Close();
		}
		
		
	}
}