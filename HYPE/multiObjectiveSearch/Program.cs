/*
 * Created by SharpDevelop.
 * User: pc
 * Date: 11/30/2017
 * Time: 11:40 PM
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Collections.Generic;
using System.IO;

namespace multiObjectiveSearch
{
	class Program
	{
		public static void Main(string[] args)
		{
			
			HypE h = new HypE("settings.txt");
			List<chromosome> ansh = h.SearchDesignSpace();
			
			
			//print answer
			StreamWriter sw = null;
			try
			{
				sw = new StreamWriter("HyPE_output.txt");
			}
			catch
			{
				if(sw != null) 
					sw.Close();
			}
			if(sw != null)
			{
				for(int i = 0; i < ansh.Count; i++)
				{
					//sw.WriteLine("ans " + i.ToString() + " : " + ans[i].PrintCMDString());
					sw.WriteLine(ansh[i].PrintRawString("\t"));
				}
			}
			sw.Close();
			
		}
		
		
	}
}