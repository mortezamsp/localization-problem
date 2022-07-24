/*
 * Created by SharpDevelop.
 * User: pc
 * Date: 12/29/2017
 * Time: 4:35 PM
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using Catfood.Shapefile;


namespace multiObjectiveSearch
{
	/// <summary>
	/// Description of Math.
	/// </summary>
	public static class Geometry
	{
		//private double DistanceOfPointToLine(PointD linestart, PointD lineend, PointD point_)
		public static double DistanceOfPointToLine(PointD linestart, PointD lineend, PointD point_, bool isSegment = false)
		{
			/*
			double d = 0;
			double dx = linestart.X - lineend.X;
			double dy = linestart.Y - lineend.Y;
			d = Math.Abs(dy*point_.X - dx*point_.Y + linestart.X*lineend.Y - lineend.X*linestart.Y) / 
					Math.Sqrt(Math.Pow(dx, 2) + Math.Pow(dy, 2));
			
			return d;
			*/

			/*
			if(isSegment)
			{
				double dot1 = dotProduct(linestart, lineend, point_);
				if(dot1 > 0)
					return DistanceOfPointToPoint(lineend, point_);
				double dot2 = dotProduct(lineend, linestart, point_);
				if(dot2 > 0)
					return DistanceOfPointToPoint(linestart, point_);
			}
			
			return Math.Abs(crossProduct(linestart, lineend, point_) / DistanceOfPointToPoint(linestart, lineend));
			*/
			
			double mindist = Math.Min(DistanceOfPointToPoint(linestart, point_), DistanceOfPointToPoint(lineend, point_));
			
			bool midx1 = (point_.X <= lineend.X && point_.X >= linestart.X);
			bool midx2 = (point_.X >= lineend.X && point_.X <= linestart.X);
			bool midy1 = (point_.Y <= lineend.Y && point_.Y >= linestart.Y);
			bool midy2 = (point_.Y >= lineend.Y && point_.Y <= linestart.Y);
			
			if(linestart.X == lineend.X)
			{
				if(midx1 || midx2)
					return point_.X - linestart.X;
				else
					return mindist;
			}
			if(linestart.Y == lineend.Y)
			{
				if(midy1 || midy2)
					return point_.Y - linestart.Y;
				else
					return mindist;
			}
			if((midx1 || midx2) && (midy1 || midy2))
			{
				double m1 = (linestart.Y - lineend.Y) / (linestart.X - lineend.X);
				double m2 = 1 / m1;
				double b1 = linestart.Y - linestart.X * m1;
				double b2 = point_.Y - point_.X * m2;
				double x = (m1 * (b2 - b1)) / (Math.Pow(m1, 2) - 1);
				double y = m1 * x + b1;
				//double y2 = m2 * x + b2;
				
				if((x >= linestart.X && x <= lineend.X) || (x <= linestart.X && x >= lineend.X))
					return DistanceOfPointToPoint(new PointD(x, y), point_);
				return mindist;
			}
			else
			{
				return mindist;
			}
			
		}
		public static double dotProduct(PointD p1, PointD p2, PointD p3)
		{
			return (p2.X - p1.X) * (p3.X - p2.X) + (p2.Y - p1.Y) * (p3.Y - p2.Y);
		}
		public static double crossProduct(PointD p1, PointD p2, PointD p3)
		{
			return (p2.X - p1.X) * (p3.Y - p1.Y) + (p2.Y - p1.Y) * (p3.X - p1.X);
		}
		public static double DistanceOfPointToPoint(PointD p1, PointD p2)
		{
			return Math.Sqrt(Math.Pow(p1.X - p2.X, 2) + Math.Pow(p1.Y - p2.Y, 2));
		}
	}
}
