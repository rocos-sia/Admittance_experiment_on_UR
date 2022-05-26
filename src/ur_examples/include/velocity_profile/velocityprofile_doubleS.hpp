/**
 * @file velocityprofile_doubleS.cpp
 * @author jc
 * @brief  用于实现double s速度规划
 * @version 0.1
 * @date 2022-03-04
 * 
 * @copyright Copyright (c) 2022
 * 
 */



#ifndef KDL_MOTION_VELOCITYPROFILE_DOUBLE_H
#define KDL_MOTION_VELOCITYPROFILE_DOUBLE_H

#include "velocityprofile.hpp"


namespace KDL {



	/**
	 * A Trapezoidal VelocityProfile implementation.
	 * @ingroup Motion
	 */
class VelocityProfile_Trap : public VelocityProfile
	{
		// For "running" a motion profile :
		double a1,a2,a3; // coef. from ^0 -> ^2 of first part
		double b1,b2,b3; // of 2nd part
		double c1,c2,c3; // of 3rd part
		double duration;
		double t1,t2;

		// specification of the motion profile :
		double maxvel;
		double maxacc;
		double startpos;
		double endpos;
	public:

		VelocityProfile_Trap(double _maxvel=0,double _maxacc=0);
		// constructs motion profile class with <maxvel> and <maxacc> as parameters of the
		// trajectory.

		virtual void SetProfile(double pos1,double pos2);

		virtual void SetProfileDuration(
			double pos1,double pos2,double newduration
		);

		/** Compute trapezoidal profile at a given fraction of max velocity
			@param pos1 Position to start from
			@param pos2 Position to end at
			@param newvelocity Fraction of max velocity to use during the
			non-ramp, flat-velocity part of the profile.
			@param KDL::epsilon <= newvelocity <= 1.0 (forcibly clamped to
			this range internally)
		*/
		virtual void SetProfileVelocity(
			double pos1,double pos2,double newvelocity
		);

        virtual void SetMax(double _maxvel,double _maxacc);
		virtual double Duration() const;
		virtual double Pos(double time) const;
		virtual double Vel(double time) const;
		virtual double Acc(double time) const;
		virtual void Write(std::ostream& os) const;
		virtual VelocityProfile* Clone() const;
		// returns copy of current VelocityProfile object. (virtual constructor)
		virtual ~VelocityProfile_Trap();
	};



}


#endif
