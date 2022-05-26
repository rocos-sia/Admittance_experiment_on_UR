
#include "velocity_profile/interpolate.h"

namespace rocos
{
    double PULSE::pos( double time )
    {
        double ret = 0;

        int n    = time / _cycle;      // 第几个周期
        double m = time - n * _cycle;  // 减去正周期整数倍的剩余时间

        switch ( _type )
        {
            case VEL:  // 速度方波，则位移三角波

                if ( m > _cycle * _ratio )  // 剩余时间超过半个周期
                {
                    ret = n * ( _amp * ( _cycle * ( 2.0 * _ratio - 1.0 ) ) ) + _amp * ( _cycle * _ratio ) -
                          _amp * ( m - _cycle * _ratio );
                }
                else
                {
                    ret = n * ( _amp * ( _cycle * ( _ratio * 2.0 - 1.0 ) ) ) + _amp * m;
                }
                return ret;
                break;

            case ACC:  // 加速度方波,速度是三角波，位移是递增
                double v0 = 0;
                for ( int i = 0; i < n; i++ )
                {
                    double t = _cycle * _ratio;
                    ret      = ret + v0 * t + 0.5 * _amp * t * t;
                    v0       = v0 + _amp * t;
                    t        = _cycle - _cycle * _ratio;
                    ret      = ret + v0 * t - 0.5 * _amp * t * t;
                    v0       = v0 - _amp * t;
                }
                if ( m > _cycle * _ratio )
                {
                    ret = ret + v0 * _cycle * _ratio + 0.5 * _amp * ( _cycle * _ratio ) * ( _cycle * _ratio );
                    v0  = v0 + _amp * ( _cycle * _ratio );
                    ret = ret + v0 * ( _cycle * ( m - _cycle * _ratio ) ) -
                          0.5 * _amp * ( _cycle * ( m - _cycle * _ratio ) ) * ( _cycle * ( m - _cycle * _ratio ) );
                    v0 = v0 - _amp * ( _cycle * ( m - _cycle * _ratio ) );
                }
                else
                {
                    ret = ret + v0 * ( _cycle * m ) + 0.5 * _amp * ( m ) * ( m );
                }
                return ret;
                break;
        }

        return ret;
    }

    double PULSE::vel( double time )
    {
        double ret = 0;

        int n    = time / _cycle;      // 第几个周期
        double m = time - n * _cycle;  // 减去正周期整数倍的剩余时间

        switch ( _type )
        {
            case VEL:
                if ( m > _cycle * _ratio )
                {
                    return -_amp;
                }
                else
                {
                    return _amp;
                }
                break;

            case ACC:
                PULSE temp( VEL, _cycle, _ratio, _amp );
                return temp.pos( time );
                return ret;
                break;
        }

        return ret;
    }

// 通用插补
#pragma region universel

    bool R_INTERP::planProfile( )
    {
        if ( _pTarget == _pStart )
        {
            printf( "please input valid target!\n" );
        }

        if ( _v_limit == 0.0 || _a_limit == 0.0 )
        {
            printf( "please input valid velocity or acceleration limit!\n" );
        }
        //TODO 改动二：new->std::make_shared
        switch ( _pType )
        {
            case trapezoid:
                if ( _trapezoid == nullptr )
                {
                    _trapezoid = std::make_shared< Trapezoid >( );
                }
                _pCurProfileType = _trapezoid;
                _trapezoid->planTrapezoidProfile( _tStart, _pStart, _pTarget, _vStart, _vTarget, _v_limit, _a_limit );
                break;

            case doubleS:
                if ( _doubleS == nullptr )
                {
                    _doubleS = std::make_shared< DoubleS >( );
                }
                _pCurProfileType = _doubleS;
                _doubleS->planDoubleSProfile( _tStart, _pStart, _pTarget, _vStart, _vTarget, _v_limit, _a_limit,
                                              _j_limit );
                break;

            case polynomial_3rd_T:
                if ( _tTerminal == 0 || _tTerminal <= _tStart )
                {
                    printf( "please input valid terminal time!\n" );
                }
                if ( _polynomial == nullptr )
                {
                    _polynomial = std::make_shared< Polynomial >( );
                }
                _pCurProfileType = _polynomial;
                _polynomial->plan3rdProfileT( _tStart, _tTerminal, _pStart, _pTarget, _vStart, _vTarget );
                break;

            case polynomial_5th_T:
                if ( _tTerminal == 0 || _tTerminal <= _tStart )
                {
                    printf( "please input valid terminal time!\n" );
                }
                if ( _polynomial == nullptr )
                {
                    _polynomial = std::make_shared< Polynomial >( );
                }
                _pCurProfileType = _polynomial;
                _polynomial->plan5thProfileT( _tStart, _tTerminal, _pStart, _pTarget, _vStart, _vTarget, _aStart,
                                              _aTarget );
                break;

            case harmonic_T:
                if ( _tTerminal == 0 || _tTerminal <= _tStart )
                {
                    printf( "please input valid terminal time!\n" );
                }
                if ( _trigonometric == nullptr )
                {
                    _trigonometric = std::make_shared< Trigonometric >( );
                }
                _pCurProfileType = _trigonometric;
                _trigonometric->planHarmonicProfileT( _pStart, _pTarget, _tStart, _tTerminal );
                break;

            case cycloidal_T:
                if ( _tTerminal == 0 || _tTerminal <= _tStart )
                {
                    printf( "please input valid terminal time!\n" );
                }
                if ( _trigonometric == nullptr )
                {
                    _trigonometric = std::make_shared< Trigonometric >( );
                }
                _pCurProfileType = _trigonometric;
                _trigonometric->planCycloidalProfileT( _pStart, _pTarget, _tStart, _tTerminal );
                break;

            case elliptic_T:
                if ( _tTerminal == 0 || _tTerminal <= _tStart )
                {
                    printf( "please input valid terminal time!\n" );
                }
                if ( _trigonometric == nullptr )
                {
                    _trigonometric = std::make_shared< Trigonometric >( );
                }
                _pCurProfileType = _trigonometric;
                _trigonometric->planEllipticProfileT( _pStart, _pTarget, _tStart, _tTerminal, 1.5 );
                break;

            case unDefine:
                _pType = doubleS;
                planProfile( );
                if ( !_pCurProfileType->isValidMovement( ) )
                {
                    _pType = trapezoid;
                    planProfile( );
                }
                break;
        }

        return isValidMovement( );
    }

    bool R_INTERP::planDoubleSPorfile( double t0, double p0, double pf,
                                       double v0, double vf, double v_limit, double a_limit, double j_limit )
    {
        setStartTime( t0 );
        setStart( p0, v0 );
        setTarget( pf, vf );
        setLimit( v_limit, a_limit, j_limit );
        if ( _doubleS == nullptr )
        {
            //TODO 改动三:new ->std::make_shared
            _doubleS = std::make_shared< DoubleS >( );
        }
        _pCurProfileType = _doubleS;
        _doubleS->planDoubleSProfile( t0, p0, pf, v0, vf, v_limit, a_limit, j_limit );
        return _doubleS->isValidMovement( );
    }

    bool R_INTERP::planTrapezoidProfile( double t0, double p0, double pf,
                                         double v0, double vf, double v_limit, double a_limit )
    {
        setStart( p0, v0 );
        setTarget( pf, vf );
        setStartTime( t0 );
        setLimit( v_limit, a_limit );
        if ( _trapezoid == nullptr )
        {  //TODO 改动四:new ->std::make_shared
            _trapezoid = std::make_shared< Trapezoid >( );
        }
        _pCurProfileType = _trapezoid;
        _trapezoid->planTrapezoidProfile( t0, p0, pf, v0, vf, v_limit, a_limit );
        return _trapezoid->isValidMovement( );
    }

    bool R_INTERP::plan3rdProfileT( double t0, double tf, double p0,
                                    double pf, double v0, double vf )
    {
        setStartTime( t0 );
        setTargetTime( tf );
        setStart( p0, v0 );
        setTarget( pf, vf );

        if ( _polynomial == nullptr )
        {  //TODO 改动5:new ->std::make_shared
            _polynomial = std::make_shared< Polynomial >( );
        }
        _pCurProfileType = _polynomial;
        _polynomial->plan3rdProfileT( t0, tf, p0, pf, v0, vf );
        return _polynomial->isValidMovement( );
    }

    bool R_INTERP::plan5thProfileT( double t0, double tf, double p0,
                                    double pf, double v0, double vf, double a0, double af )
    {
        setStartTime( t0 );
        setTargetTime( tf );
        setStart( p0, v0, a0 );
        setTarget( pf, vf, af );

        if ( _polynomial == nullptr )
        {  //TODO 改动6:new ->std::make_shared
            _polynomial = std::make_shared< Polynomial >( );
        }
        _pCurProfileType = _polynomial;
        _polynomial->plan5thProfileT( t0, tf, p0, pf, v0, vf, a0, af );
        return _polynomial->isValidMovement( );
    }

#pragma endregion universel

// T型波
#pragma region TRAPEZOID

    void Trapezoid::planProfile( double t0, double p0, double pf, double v0, double vf, double v_limit, double a_limit,
                                 double j_limit )
    {
        planTrapezoidProfile( t0, p0, pf, v0, vf, v_limit, a_limit );
    }

    void Trapezoid::planTrapezoidProfile( double t0, double p0, double pf,
                                          double v0, double vf, double v_limit, double a_limit )
    {
        _plannedProfile = false;

        _direction = sign( pf - p0 );

        if ( _direction == 0 )
        {
            //		std::cerr<<("target and start points are the same!\n");
            return;
        }

        double p0_z = _direction * p0;
        double pf_z = _direction * pf;
        double v0_z = _direction * v0;
        double vf_z = _direction * vf;

        double vmax = fabs( v_limit );
        double amax = fabs( a_limit );
        double vmin = -vmax;
        double amin = -amax;

        double v_limit_z = ( _direction + 1.0 ) / 2.0 * vmax + ( _direction - 1.0 ) / 2.0 * vmin;
        double a_limit_z = ( _direction + 1.0 ) / 2.0 * amax + ( _direction - 1.0 ) / 2.0 * amin;

        trapezoid_check_T( t0, p0_z, pf_z, v0_z, vf_z, v_limit_z, a_limit_z );
    }

    void Trapezoid::trapezoid_check_T( double t0, double p0, double pf,
                                       double v0, double vf, double v_limit, double a_limit )
    {
        // the 'trapezoid' profile trajectory looks like :
        //
        //	   v |	_________
        //		 | /         \
	    //       |/           \
	    //       |             \
	    //	     |              |
        //       0-1--------2---3--> t
        // let Ta Tv Td be the acceleration cursing and deceleration time interval
        // so, pf - p0 = v0*Ta +0.5*_a[1]*Ta*Ta + _v[1]*Tv + vf*Td - 0.5*_a[3]*Td*Td
        // Ta = (_v[1] - _v[0])/_a[1];	// here, _a[1],acceleration, bigger than zero
        // Td = (_v[3] - _v[2])/_a[3];	// here, _a[3],deceleration, samller than zero,
        // Tv = ((pf-p0) - (_v[0]*Ta + 0.5*_a[1]*Ta*Ta) - (_v[3]*Td - 0.5*_a[3]*Td*Td))/_v[1];

        _t[ 0 ] = t0;  // start time
        _x[ 0 ] = p0;
        _x[ 3 ] = pf;
        _v[ 0 ] = v0;
        _v[ 3 ] = vf;
        _a[ 0 ] = 0;
        _a[ 2 ] = 0;

        double amax = a_limit;  // a_limit: acceleration in phase #1;
        _a[ 1 ]     = a_limit;  // plan the trajectory with the limited acceleration

        double amin = -a_limit;  //amin: acceleration in phase #3, be deceleration.
        _a[ 3 ]     = amin;      // by default, deceleration equals to minus acceleration

        // it is necessary to check whether the trajectory is feasible or not by
        if ( a_limit * ( pf - p0 ) < fabs( v0 * v0 - vf * vf ) / 2.0 )
        {
            std::cerr << ( "trapezoid trajectory is not feasible\n" );
            _plannedProfile = false;
            return;
        }
        else
        {
            double Ta, Tv, Td;
            if ( ( pf - p0 ) * a_limit > v_limit * v_limit - ( v0 * v0 + vf * vf ) / 2.0 )
            {
                // Case 1: _v[1] = v_limit
                _v[ 1 ] = v_limit;
            }
            else
            {
                // Case 2: _v[1] < v_limit
                _v[ 1 ] = sqrt( ( pf - p0 ) * a_limit + ( v0 * v0 + vf * vf ) / 2.0 );
            }

            _v[ 2 ] = _v[ 1 ];

            Ta = ( _v[ 1 ] - _v[ 0 ] ) / _a[ 1 ];
            Td = ( _v[ 3 ] - _v[ 2 ] ) / _a[ 3 ];
            Tv = ( ( pf - p0 ) - ( _v[ 0 ] * Ta + 0.5 * _a[ 1 ] * Ta * Ta ) - ( _v[ 3 ] * Td - 0.5 * _a[ 3 ] * Td * Td ) ) / _v[ 1 ];

            _t[ 1 ] = _t[ 0 ] + Ta;
            _t[ 2 ] = _t[ 1 ] + Tv;
            _t[ 3 ] = _t[ 2 ] + Td;

            _x[ 1 ] = _x[ 0 ] + _v[ 0 ] * Ta / 2.0 + _v[ 1 ] * ( Ta - Ta / 2.0 );
            _x[ 2 ] = _x[ 0 ] + _v[ 0 ] * Ta / 2.0 + _v[ 1 ] * ( Tv - Ta / 2.0 );
        }

        _plannedProfile = true;
    }

    double Trapezoid::scaleToDuration( double newDuration )
    {
        if ( newDuration < ( _t[ 3 ] - _t[ 0 ] ) )
        {
            //            std::cerr<<("New duration must be longer!\n");
            return _t[ 3 ] - _t[ 0 ];
        }

        double p0 = _x[ 0 ];
        double pf = _x[ 3 ];
        double v0 = _v[ 0 ];
        double vf = _v[ 3 ];
        double h  = pf - p0;

        double T  = _t[ 3 ] - _t[ 0 ];
        double Ta = _t[ 1 ] - _t[ 0 ];
        double Tv = _t[ 2 ] - _t[ 1 ];
        double Td = _t[ 3 ] - _t[ 2 ];

        Ta = Ta * newDuration / T;
        Tv = Tv * newDuration / T;
        Td = Td * newDuration / T;

        // 3 equations to get new _v[1] _a[1] and _a[3]
        // Ta = (_v[1] - v0)/_a[1];
        // Tv = ((pf-p0) - (_v[0]*Ta + 0.5*_a[1]*Ta*Ta) - (_v[3]*Td - 0.5*_a[3]*Td*Td))/_v[1];
        // Td = (_v[3] - _v[1])/_a[3];

        _a[ 1 ] = -( 1.0 * ( 2.0 * p0 - 2.0 * pf + 2.0 * Ta * v0 + Td * v0 + 2.0 * Tv * v0 + Td * vf ) ) /
                  ( Ta * ( Ta + Td + 2.0 * Tv ) );
        _a[ 3 ] = ( 2.0 * p0 - 2.0 * pf + Ta * v0 + Ta * vf + 2.0 * Td * vf + 2.0 * Tv * vf ) / ( Td * ( Ta + Td + 2.0 * Tv ) );
        _v[ 1 ] = -( 1.0 * ( 2.0 * p0 - 2.0 * pf + Ta * v0 + Td * vf ) ) / ( Ta + Td + 2.0 * Tv );

        _v[ 2 ] = _v[ 1 ];

        Ta = ( _v[ 1 ] - _v[ 0 ] ) / _a[ 1 ];
        Td = ( _v[ 3 ] - _v[ 2 ] ) / _a[ 3 ];
        Tv = ( ( pf - p0 ) - ( _v[ 0 ] * Ta + 0.5 * _a[ 1 ] * Ta * Ta ) - ( _v[ 3 ] * Td - 0.5 * _a[ 3 ] * Td * Td ) ) / _v[ 1 ];

        _t[ 1 ] = _t[ 0 ] + Ta;
        _t[ 2 ] = _t[ 1 ] + Tv;
        _t[ 3 ] = _t[ 2 ] + Td;

        _x[ 1 ] = _x[ 0 ] + _v[ 0 ] * Ta / 2.0 + _v[ 1 ] * ( Ta - Ta / 2.0 );
        _x[ 2 ] = _x[ 0 ] + _v[ 0 ] * Ta / 2.0 + _v[ 1 ] * ( Tv - Ta / 2.0 );

        return getDuration( );
    }

    bool Trapezoid::isValidMovement( ) const
    {
        return _plannedProfile;
    }

    double Trapezoid::getDuration( ) const
    {
        return ( _t[ 3 ] - _t[ 0 ] );
    }

    double Trapezoid::pos( double t ) const
    {
        double p;

        double dt = ( t - _t[ 0 ] );

        double ta = _t[ 1 ] - _t[ 0 ];
        double tv = _t[ 2 ] - _t[ 1 ];
        double td = _t[ 3 ] - _t[ 2 ];
        double T  = _t[ 3 ] - _t[ 0 ];

        if ( dt < 0 )
        {
            p = _x[ 0 ];
        }

        if ( dt < ta && dt >= 0 )
        {
            p = _x[ 0 ] + _v[ 0 ] * dt + ( _v[ 1 ] - _v[ 0 ] ) / 2.0 / ta * dt * dt;
        }

        if ( ta <= dt && dt < ta + tv )
        {
            p = _x[ 0 ] + _v[ 0 ] * ta / 2.0 + _v[ 1 ] * ( dt - ta / 2.0 );
        }

        if ( ta + tv <= dt && dt < T )
        {
            p = _x[ 3 ] - _v[ 3 ] * ( T - dt ) - ( _v[ 2 ] - _v[ 3 ] ) / 2.0 / td * ( T - dt ) * ( T - dt );
        }

        if ( dt >= T )
        {
            dt = T;
            p  = _x[ 3 ] - _v[ 3 ] * ( T - dt ) - ( _v[ 2 ] - _v[ 3 ] ) / 2.0 / td * ( T - dt ) * ( T - dt );
        }

        p = p * _direction;

        return p;
    }

    double Trapezoid::vel( double t ) const
    {
        double v;

        double dt = ( t - _t[ 0 ] );

        double ta = _t[ 1 ] - _t[ 0 ];
        double tv = _t[ 2 ] - _t[ 1 ];
        double td = _t[ 3 ] - _t[ 2 ];
        double T  = _t[ 3 ] - _t[ 0 ];

        if ( dt < 0 )
        {
            v = _v[ 0 ];
        }

        if ( dt < ta && dt >= 0 )
        {
            v = _v[ 0 ] + ( _v[ 1 ] - _v[ 0 ] ) / ta * dt;
        }

        if ( ta <= dt && dt < ta + tv )
        {
            v = _v[ 1 ];
        }

        if ( ta + tv <= dt && dt < T )
        {
            v = _v[ 3 ] + ( _v[ 2 ] - _v[ 3 ] ) / td * ( T - dt );
        }

        if ( dt >= T )
        {
            dt = T;
            v  = _v[ 3 ] + ( _v[ 2 ] - _v[ 3 ] ) / td * ( T - dt );
        }

        v = v * _direction;

        return v;
    }

    double Trapezoid::acc( double t ) const
    {
        double a;
        double dt = ( t - _t[ 0 ] );

        double ta = _t[ 1 ] - _t[ 0 ];
        double tv = _t[ 2 ] - _t[ 1 ];
        double td = _t[ 3 ] - _t[ 2 ];
        double T  = _t[ 3 ] - _t[ 0 ];

        if ( dt < 0 )
        {
            a = 0.0;
        }

        if ( dt < ta && dt >= 0 )
        {
            a = _a[ 1 ];
        }

        if ( ta <= dt && dt < ta + tv )
        {
            a = _a[ 2 ];
        }

        if ( ta + tv <= dt && dt < T )
        {
            a = _a[ 3 ];
        }

        if ( dt >= T )
        {
            a = 0.0;
        }

        a = a * _direction;
        return a;
    }

    double Trapezoid::jerk( double t ) const
    {
        return 0;
    }

    double Trapezoid::max_vel( ) const
    {
        double v = 0;

        for ( int i = 0; i < 4; i++ )
        {
            v = fmax( fabs( v ), fabs( _v[ i ] ) );
        }

        for ( int i = 0; i < 4; i++ )
        {
            if ( v == fabs( _v[ i ] ) )
            {
                v = _v[ i ] * _direction;
            }
        }
        return v;
    }

    double Trapezoid::max_acc( ) const
    {
        double a = 0;
        for ( int i = 0; i < 4; i++ )
        {
            a = fmax( fabs( a ), fabs( _a[ i ] ) );
        }
        for ( int i = 0; i < 4; i++ )
        {
            if ( fabs( a ) == fabs( _a[ i ] ) )
            {
                a = _a[ i ] * _direction;
            }
        }
        return a;
    }

#pragma endregion TRAPEZOID

// S型波
#pragma region DOUBLES

    void DoubleS::planProfile( double t0, double p0, double pf, double v0, double vf, double v_limit, double a_limit,
                               double j_limit )
    {
        planDoubleSProfile( t0, p0, pf, v0, vf, v_limit, a_limit, j_limit );
    }

    void DoubleS::planDoubleSProfile( double t0, double p0, double pf,
                                      double v0, double vf, double v_limit, double a_limit, double j_limit )
    {
        _plannedProfile = false;

        _direction = sign( pf - p0 );

        if ( _direction == 0 )
        {
            _t[ 0 ] = t0;
            _t[ 7 ] = t0;
            std::cerr << ( "target and start points are the same!\n" );
            return;
        }

        double p0_z = _direction * p0;
        double pf_z = _direction * pf;
        double v0_z = _direction * v0;
        double vf_z = _direction * vf;

        double vmax = fabs( v_limit );
        double amax = fabs( a_limit );
        double jmax = fabs( j_limit );

        double vmin = -vmax;
        double amin = -amax;
        double jmin = -jmax;

        double v_limit_z = ( _direction + 1.0 ) / 2.0 * vmax + ( _direction - 1.0 ) / 2.0 * vmin;
        double a_limit_z = ( _direction + 1.0 ) / 2.0 * amax + ( _direction - 1.0 ) / 2.0 * amin;
        double j_limit_z = ( _direction + 1.0 ) / 2.0 * jmax + ( _direction - 1.0 ) / 2.0 * jmin;

        doubleS_check_T( t0, p0_z, pf_z, v0_z, vf_z, v_limit_z, a_limit_z, j_limit );
    }

    void DoubleS::doubleS_check_T( double t0, double p0, double pf,
                                   double v0, double vf, double v_limit, double a_limit, double j_limit )
    {
        // the 'doubleS' has acceleration phase, cursing phase and deceleration phase
        // in acceleration phase and deceleration phase, the acceleration has three phases:
        //			acc acceleration phase, acc cursing phase and acc deceleration phase
        // a |	_________
        //   | /         \
	//       |/           \
	//       0--1- -----2--3-------4--5------6--7-----------> t
        //								\			/
        //								 \		   /
        //								  ---------
        // the time interval between 0~1 is: Ta_acc
        // the time interval between 1~2 is: Tv_acc
        // the time interval between 2~3 is: Td_acc
        // the time interval between 0~3 is: Ta

        // the time interval between 3~4 is: Tv

        // the time interval between 4~5 is: Ta_dec
        // the time interval between 5~6 is: Tv_dec
        // the time interval between 6~7 is: Td_dec
        // the time interval between 4~7 is: Td

        // by default, Td_acc = Ta_acc, and Ta_dec = Td_dec

        _t[ 0 ] = t0;  // start time

        _x[ 0 ] = p0;
        _x[ 7 ] = pf;

        _v[ 0 ] = v0;
        _v[ 7 ] = vf;

        _a[ 0 ] = 0;
        _a[ 4 ] = 0;
        _a[ 3 ] = 0;
        _a[ 7 ] = 0;

        _j[ 0 ] = 0;
        _j[ 2 ] = 0;
        _j[ 4 ] = 0;
        _j[ 6 ] = 0;

        double r = 1;  // scale factor for re-planning doubleS profile

        // It is necessary to check if it is possible to perform the trajectory
        // with a double jerk impulse(one positive and one negative) only.
        // For this purpose, let us define
        double Tjx = fmin( sqrt( fabs( vf - v0 ) / j_limit ), a_limit / j_limit );
        if ( Tjx == a_limit / j_limit )
        {
            if ( pf - p0 <= 0.5 * ( v0 + vf ) * ( Tjx + fabs( vf - v0 ) / a_limit ) )
            {
                _plannedProfile = false;
                std::cerr << ( "The two points are two close! DoubleS trajectory is not feasible!\n" );
                return;
            }
        }
        else
        {
            if ( pf - p0 <= Tjx * ( vf + v0 ) )
            {
                _plannedProfile = false;
                std::cerr << ( "The two points are two close! DoubleS trajectory is not feasible!\n" );
                return;
            }
        }

        // In this case, by defining the maximum value of the velocity during
        // the motion as vmax = max(v(t)), there are two possibilities:
        // Case 1: vmax = v_limit
        // Case 2: vmax < v_limit
        // In both Case 1 and Case 2, it is possible that the limited
        // acceleration(positive or negative) is not reached. This may happen if
        // the displacement is small, if the limited allowed acceleration a_limit is high.

        // Let us define:
        // Ta:		acceleration period;
        // Tv:		constant velocity period;
        // Td:		deceleration period;
        // T:		total duration of the trajectory (=Ta+Tv+Td);
        // Ta_acc:	acceleration increasing time-interval during acceleration peroid;
        // Tv_acc:	acceleration keep constant time-interval during acceleration period;
        // Td_acc:  acceleration decreasing time-interval during acceleration peroid;
        // Ta_dec:  deceleration increasing time-interval during deceleration peroid;
        // Tv_dec:  deceleration keep constant time-interval during decelration peroid;
        // Td_dec:  deceleration decreasing time-interval during deceleration period.

        double Ta, Tv, Td, T, Ta_acc, Tv_acc, Td_acc, Ta_dec, Tv_dec, Td_dec;
        double jmax_a, jmin_a, jmax_d, jmin_d, amax_a, amin_d, vmax;

        // Case 1: vmax = v_limit
        if ( ( v_limit - v0 ) * j_limit < a_limit * a_limit )  // acceleration phase
        {
            // a_limit is not reached
            Ta_acc = sqrt( ( v_limit - v0 ) / j_limit );
            Ta     = 2 * Ta_acc;
        }
        else
        {
            // a_limit is reached
            Ta_acc = a_limit / j_limit;
            Ta     = Ta_acc + ( v_limit - v0 ) / a_limit;
        }

        if ( ( v_limit - vf ) * j_limit < a_limit * a_limit )  // deceleration phase
        {
            // a_limit is not reached
            Ta_dec = sqrt( ( v_limit - vf ) / j_limit );
            Td     = 2 * Ta_dec;
        }
        else
        {
            // a_limit is reached
            Ta_dec = a_limit / j_limit;
            Td     = Ta_dec + ( v_limit - vf ) / a_limit;
        }

        // Finally, it is possible to determine the time duration of the
        // constant velocity segments as
        Tv = ( pf - p0 ) / v_limit - Ta / 2.0 * ( 1.0 + v0 / v_limit ) - Td / 2.0 * ( 1.0 + vf / v_limit );

        if ( Tv <= 0 )
        {
            // Case 2: vmax < v_max
            // In this case, the constant velocity segment is not present(Tv=0),
            // and the duration of the acceleration and deceleration segments
            // can be easily computed if the limited accelerations are reached
            // in both segments. In this case
            double temp = 0;
            while ( 1 )
            {
                Ta_acc = a_limit / j_limit;
                Ta_dec = a_limit / j_limit;
                temp   = ( a_limit * a_limit * a_limit * a_limit ) / ( j_limit * j_limit ) + 2.0 * ( v0 * v0 + vf * vf ) + a_limit * ( 4.0 * ( pf - p0 ) - 2.0 * a_limit / j_limit * ( v0 + vf ) );
                Ta     = ( ( a_limit * a_limit ) / j_limit - 2.0 * v0 + sqrt( temp ) ) / 2.0 / a_limit;
                Td     = ( ( a_limit * a_limit ) / j_limit - 2.0 * vf + sqrt( temp ) ) / 2.0 / a_limit;
                Tv     = 0;
                // if Ta<2*Ta_acc or Td<2*Ta_dec, then the limited acceleration is
                // not reached in at least one of the two segments, and therefore it
                // is not possible to use above eq. In this case(indeed rather
                // unusual), the determination of parameters is quite difficult, and
                // it may be more convenient to find an approximated solution that,
                // although not optimal, results acceptable from a computational
                // point of view. A possible way to determine this solution is to
                // progressively decrease the value of a_limit (e.g. by assuming
                // a_limit=r*a_limit, with 0<r<1) and compute the duration of the
                // segments by means of above eq. until the condition Ta>2*Ta_acc and
                // Td>2*Ta_dec are both true.

                if ( Ta < 0 || Td < 0 )
                {
                    // During this recursive computation,it may happen that Ta
                    // or Td became negative. In this case, only one of the
                    // acceleration or deceleration phase is necessary,
                    // depending on the values of the initial and final
                    // velocities. If Ta<0 (note that in this case necessarily
                    // v0>v1), the acceleration phase is not present. Then, Ta
                    // is set to 0 and the duration of the deceleration segment
                    // can be computed according to

                    if ( Ta < 0 )
                    {
                        Td     = 2 * ( pf - p0 ) / ( vf + v0 );
                        Ta_dec = ( j_limit * ( pf - p0 ) - sqrt( j_limit * ( j_limit * ( pf - p0 ) * ( pf - p0 ) +
                                                                             ( vf + v0 ) * ( vf + v0 ) * ( vf - v0 ) ) ) ) /
                                 ( j_limit * ( vf + v0 ) );
                        Ta     = 0.0;
                        Ta_acc = 0.0;
                        break;
                    }

                    // In the dual case, i.e. when Td<0 (this case is possible
                    // when v1>v0), the deceleration phase is not necessary
                    // (Td=0) and the duration of the acceleration period must
                    // be computed as
                    if ( Td < 0 )
                    {
                        Ta     = 2 * ( pf - p0 ) / ( vf + v0 );
                        Ta_acc = ( j_limit * ( pf - p0 ) - sqrt( j_limit * ( j_limit * ( pf - p0 ) * ( pf - p0 ) -
                                                                             ( vf + v0 ) * ( vf + v0 ) * ( vf - v0 ) ) ) ) /
                                 ( j_limit * ( vf + v0 ) );
                        Td     = 0;
                        Ta_dec = 0;
                        break;
                    }
                }
                else
                {
                    if ( ( Ta < 2 * Ta_acc ) || ( Td < 2 * Ta_dec ) )
                    {
                        r       = r - 0.01;
                        a_limit = r * a_limit;
                        //std::cerr<<("DoubleS re-planning: Ta<2*Tj1 or Td<2*Tj2\n");
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        amax_a = j_limit * Ta_acc;
        amin_d = -j_limit * Ta_dec;
        vmax   = v0 + ( Ta - Ta_acc ) * amax_a;

        jmax_a = j_limit;
        jmin_a = -j_limit;
        jmax_d = -j_limit;
        jmin_d = j_limit;

        Td_acc = Ta_acc;
        Tv_acc = Ta - 2.0 * Ta_acc;

        Td_dec = Ta_dec;
        Tv_dec = Td - 2.0 * Ta_dec;

        T = Ta + Tv + Td;

        _j[ 1 ] = jmax_a;
        _j[ 3 ] = jmin_a;
        _j[ 5 ] = jmax_d;
        _j[ 7 ] = jmin_d;

        _t[ 0 ] = t0;
        _t[ 1 ] = _t[ 0 ] + Ta_acc;
        _t[ 2 ] = _t[ 1 ] + Tv_acc;
        _t[ 3 ] = _t[ 2 ] + Td_acc;
        _t[ 4 ] = _t[ 3 ] + Tv;
        _t[ 5 ] = _t[ 4 ] + Ta_dec;
        _t[ 6 ] = _t[ 5 ] + Tv_dec;
        _t[ 7 ] = _t[ 6 ] + Td_dec;

        _a[ 1 ] = amax_a;
        _a[ 2 ] = amax_a;
        _a[ 5 ] = amin_d;
        _a[ 6 ] = amin_d;

        _v[ 3 ] = vmax;
        _v[ 4 ] = vmax;

        _v[ 1 ] = vel( _t[ 1 ] );
        _v[ 2 ] = vel( _t[ 2 ] );
        _v[ 5 ] = vel( _t[ 5 ] );
        _v[ 6 ] = vel( _t[ 6 ] );

        _plannedProfile = true;
    }

    double DoubleS::scaleToDuration( double newDuration )
    {
        if ( newDuration <= ( _t[ 7 ] - _t[ 0 ] ) )
        {
            std::cerr << ( "New duration must be longer!\n" );
            return _t[ 7 ] - _t[ 0 ];
        }

        double t0 = _t[ 0 ];
        double p0 = _x[ 0 ];
        double pf = _x[ 7 ];
        double v0 = _v[ 0 ];
        double vf = _v[ 7 ];

        double j_limit = _j[ 1 ];

        double Ta = _t[ 3 ] - _t[ 0 ];
        double Tv = _t[ 4 ] - _t[ 3 ];
        double Td = _t[ 7 ] - _t[ 4 ];
        double T  = _t[ 7 ] - _t[ 0 ];

        Ta = Ta * newDuration / T;
        Tv = Tv * newDuration / T;
        Td = Td * newDuration / T;

        // decrease the cursing velocity
        double vmax   = -( p0 - pf + ( Ta * v0 ) / 2 + ( Td * vf ) / 2 ) / ( Ta / 2 + Td / 2 + Tv );
        double amax_a = _a[ 1 ];

        // decrease jmax_a jmin_a jmax_d jmin_d and re-plan time interval
        double Ta_acc, Ta_dec, jmax_a, jmax_d, amin_d;
        if ( ( vmax - v0 ) * j_limit < amax_a * amax_a )
        {
            Ta_acc = 0.5 * Ta;
            jmax_a = ( vmax - v0 ) / ( Ta_acc * Ta_acc );
        }
        else
        {
            Ta_acc = Ta - ( vmax - v0 ) / amax_a;
            jmax_a = amax_a / Ta_acc;
        }

        if ( ( vmax - vf ) * j_limit < amax_a * amax_a )
        {
            Ta_dec = 0.5 * Td;
            jmax_d = -( vmax - vf ) / ( Ta_dec * Ta_dec );
        }
        else
        {
            Ta_dec = Td - ( vmax - vf ) / amax_a;
            jmax_d = -amax_a / Ta_dec;
        }

        amax_a = jmax_a * Ta_acc;
        amin_d = jmax_d * Ta_dec;
        vmax   = v0 + ( Ta - Ta_acc ) * amax_a;

        double tempTv = ( pf - p0 ) / vmax - Ta / 2.0 * ( 1.0 + v0 / vmax ) - Td / 2.0 * ( 1.0 + vf / vmax );
        if ( tempTv < 0 || fabs( tempTv - Tv ) > 1e-4 )
        {
            std::cerr << ( "Re-plan to scale duration failed!\n" );
            return getDuration( );
        }

        T = Ta + Tv + Td;

        _j[ 1 ] = jmax_a;
        _j[ 3 ] = -jmax_a;
        _j[ 5 ] = jmax_d;
        _j[ 7 ] = -jmax_d;

        double Td_acc = Ta_acc;
        double Tv_acc = Ta - 2.0 * Ta_acc;

        double Td_dec = Ta_dec;
        double Tv_dec = Td - 2.0 * Ta_dec;

        _t[ 0 ] = t0;
        _t[ 1 ] = _t[ 0 ] + Ta_acc;
        _t[ 2 ] = _t[ 1 ] + Tv_acc;
        _t[ 3 ] = _t[ 2 ] + Td_acc;
        _t[ 4 ] = _t[ 3 ] + Tv;
        _t[ 5 ] = _t[ 4 ] + Ta_dec;
        _t[ 6 ] = _t[ 5 ] + Tv_dec;
        _t[ 7 ] = _t[ 6 ] + Td_dec;

        _a[ 1 ] = amax_a;
        _a[ 2 ] = amax_a;
        _a[ 5 ] = amin_d;
        _a[ 6 ] = amin_d;

        _v[ 3 ] = vmax;
        _v[ 4 ] = vmax;

        _v[ 1 ] = vel( _t[ 1 ] );
        _v[ 2 ] = vel( _t[ 2 ] );
        _v[ 5 ] = vel( _t[ 5 ] );
        _v[ 6 ] = vel( _t[ 6 ] );

        _plannedProfile = true;

        return getDuration( );
    }

    double DoubleS::scaleToDuration_unlimited( double newDuration )
    {
        double t0 = _t[ 0 ];
        double p0 = _x[ 0 ];
        double pf = _x[ 7 ];
        double v0 = _v[ 0 ];
        double vf = _v[ 7 ];

        double j_limit = _j[ 1 ];

        double Ta = _t[ 3 ] - _t[ 0 ];
        double Tv = _t[ 4 ] - _t[ 3 ];
        double Td = _t[ 7 ] - _t[ 4 ];
        double T  = _t[ 7 ] - _t[ 0 ];

        Ta = Ta * newDuration / T;
        Tv = Tv * newDuration / T;
        Td = Td * newDuration / T;

        // decrease the cursing velocity
        double vmax   = -( p0 - pf + ( Ta * v0 ) / 2 + ( Td * vf ) / 2 ) / ( Ta / 2 + Td / 2 + Tv );
        double amax_a = _a[ 1 ];

        // decrease jmax_a jmin_a jmax_d jmin_d and re-plan time interval
        double Ta_acc, Ta_dec, jmax_a, jmax_d, amin_d;
        if ( ( vmax - v0 ) * j_limit < amax_a * amax_a )
        {
            Ta_acc = 0.5 * Ta;
            jmax_a = ( vmax - v0 ) / ( Ta_acc * Ta_acc );
        }
        else
        {
            Ta_acc = Ta - ( vmax - v0 ) / amax_a;
            jmax_a = amax_a / Ta_acc;
        }

        if ( ( vmax - vf ) * j_limit < amax_a * amax_a )
        {
            Ta_dec = 0.5 * Td;
            jmax_d = -( vmax - vf ) / ( Ta_dec * Ta_dec );
        }
        else
        {
            Ta_dec = Td - ( vmax - vf ) / amax_a;
            jmax_d = -amax_a / Ta_dec;
        }

        amax_a = jmax_a * Ta_acc;
        amin_d = jmax_d * Ta_dec;
        vmax   = v0 + ( Ta - Ta_acc ) * amax_a;

        double tempTv = ( pf - p0 ) / vmax - Ta / 2.0 * ( 1.0 + v0 / vmax ) - Td / 2.0 * ( 1.0 + vf / vmax );
        if ( tempTv < 0 || fabs( tempTv - Tv ) > 1e-4 )
        {
            std::cerr << ( "Re-plan to scale duration failed!\n" );
            return getDuration( );
        }

        T = Ta + Tv + Td;

        _j[ 1 ] = jmax_a;
        _j[ 3 ] = -jmax_a;
        _j[ 5 ] = jmax_d;
        _j[ 7 ] = -jmax_d;

        double Td_acc = Ta_acc;
        double Tv_acc = Ta - 2.0 * Ta_acc;

        double Td_dec = Ta_dec;
        double Tv_dec = Td - 2.0 * Ta_dec;

        _t[ 0 ] = t0;
        _t[ 1 ] = _t[ 0 ] + Ta_acc;
        _t[ 2 ] = _t[ 1 ] + Tv_acc;
        _t[ 3 ] = _t[ 2 ] + Td_acc;
        _t[ 4 ] = _t[ 3 ] + Tv;
        _t[ 5 ] = _t[ 4 ] + Ta_dec;
        _t[ 6 ] = _t[ 5 ] + Tv_dec;
        _t[ 7 ] = _t[ 6 ] + Td_dec;

        _a[ 1 ] = amax_a;
        _a[ 2 ] = amax_a;
        _a[ 5 ] = amin_d;
        _a[ 6 ] = amin_d;

        _v[ 3 ] = vmax;
        _v[ 4 ] = vmax;

        _v[ 1 ] = vel( _t[ 1 ] );
        _v[ 2 ] = vel( _t[ 2 ] );
        _v[ 5 ] = vel( _t[ 5 ] );
        _v[ 6 ] = vel( _t[ 6 ] );

        _plannedProfile = true;

        return getDuration( );
    }

    double DoubleS::my_scaleToDuration( double newDuration, double t0, double p0, double pf,
                                        double v0, double vf )
    {
        if ( newDuration == 0 ) return getDuration( );
        double k = getDuration( ) / newDuration;

        planDoubleSProfile( t0, p0, pf, v0, vf, k * max_vel( ), k * k * max_acc( ), k * k * k * max_jerk( ) );
        return newDuration;
    }

    bool DoubleS::isValidMovement( ) const
    {
        return _plannedProfile;
    }

    double DoubleS::getDuration( ) const
    {
        return ( _t[ 7 ] - _t[ 0 ] );
    }

    double DoubleS::pos( double t ) const
    {
        double p;

        double dt = ( t - _t[ 0 ] );

        double Tj1, Tj2, Ta, Tv, Td, T, v_lim, a_lima, a_limd, j_limit, jmin, p0, pf, v0, vf;

        j_limit = _j[ 1 ];
        jmin    = -_j[ 7 ];

        a_lima = _a[ 1 ];
        a_limd = _a[ 5 ];

        v_lim = _v[ 3 ];

        Tj1 = _t[ 1 ] - _t[ 0 ];
        Ta  = _t[ 3 ] - _t[ 0 ];
        Tv  = _t[ 4 ] - _t[ 3 ];
        T   = _t[ 7 ] - _t[ 0 ];
        Td  = T - Tv - Ta;
        Tj2 = a_limd / jmin;

        p0 = _x[ 0 ];
        pf = _x[ 7 ];
        v0 = _v[ 0 ];
        vf = _v[ 7 ];

        if ( dt <= 0 )
        {
            p = p0;
        }

        // a) 加速度增大
        if ( dt > 0 && dt < Tj1 )
        {
            p = p0 + v0 * dt + j_limit * dt * dt * dt / 6.0;
        }

        // b) 匀加速度
        if ( dt >= Tj1 && dt < Ta - Tj1 )
        {
            p = p0 + v0 * dt + a_lima * ( 3.0 * dt * dt - 3.0 * Tj1 * dt + Tj1 * Tj1 ) / 6.0;
        }
        // c) 加速度减小
        if ( dt >= Ta - Tj1 && dt < Ta )
        {
            p = p0 + ( v_lim + v0 ) * Ta / 2 - v_lim * ( Ta - dt ) + j_limit * ( Ta - dt ) * ( Ta - dt ) * ( Ta - dt ) / 6.0;
        }
        // 2.匀速段
        // a) 匀速度
        if ( dt >= Ta && dt <= Ta + Tv )
        {
            p = p0 + ( v_lim + v0 ) * Ta / 2.0 + v_lim * ( dt - Ta );
        }
        // 3.减速段
        // a) 加速度减小
        if ( dt > T - Td && dt <= T - Td + Tj2 )
        {
            p = pf - ( v_lim + vf ) * Td / 2.0 + v_lim * ( dt - T + Td ) +
                jmin * ( dt - T + Td ) * ( dt - T + Td ) * ( dt - T + Td ) / 6.0;
            //dt = dt - Ta - Tv;
            //p = pos(Ta + Tv) + v_lim*dt + jmin*dt*dt*dt/6.0;
        }
        // b) 匀加速度
        if ( dt > T - Td + Tj2 && dt <= T - Tj2 )
        {
            p = pf - ( v_lim + vf ) * Td / 2.0 + v_lim * ( dt - T + Td ) +
                a_limd / 6.0 * ( 3.0 * ( dt - T + Td ) * ( dt - T + Td ) - 3.0 * Tj2 * ( dt - T + Td ) + Tj2 * Tj2 );
            //dt = dt - Ta - Tv - Tj2;
            //p = pos(T-Td+Tj2) + vel(Ta + Tv + Tj2)*dt + 0.5*a_limd*(dt*dt);
        }
        // c) 加速度增大
        if ( dt > T - Tj2 )
        {
            p = pf - vf * ( T - dt ) + jmin / 6.0 * ( T - dt ) * ( T - dt ) * ( T - dt );
            //dt = dt - T + Tj2;
            //p = pos(T-Tj2) + vel(T-Tj2)*dt + 0.5*acc(T-Tj2)*(dt*dt) - (jmin*dt*dt*dt)/6.0;
        }

        if ( dt >= T )
        {
            p = pf;
        }

        p = p * _direction;

        return p;
    }

    double DoubleS::vel( double t ) const
    {
        double v;

        double dt = ( t - _t[ 0 ] );

        double Tj1, Tj2, Ta, Tv, Td, T, v_lim, a_lima, a_limd, j_limit, jmin, p0, pf, v0, vf;

        j_limit = _j[ 1 ];
        //jmin = -j_limit;
        jmin = -_j[ 7 ];

        a_lima = _a[ 1 ];
        a_limd = _a[ 5 ];

        v_lim = _v[ 3 ];

        Tj1 = _t[ 1 ] - _t[ 0 ];
        Ta  = _t[ 3 ] - _t[ 0 ];
        Tv  = _t[ 4 ] - _t[ 3 ];
        T   = _t[ 7 ] - _t[ 0 ];
        Td  = T - Tv - Ta;
        Tj2 = a_limd / jmin;

        p0 = _x[ 0 ];
        pf = _x[ 7 ];
        v0 = _v[ 0 ];
        vf = _v[ 7 ];

        // a) 加速度增大
        if ( dt < Tj1 )
        {
            v = v0 + j_limit * dt * dt / 2.0;
        }

        // b) 匀加速度
        if ( dt >= Tj1 && dt < Ta - Tj1 )
        {
            v = v0 + a_lima * ( dt - Tj1 / 2 );
        }
        // c) 加速度减小
        if ( dt >= Ta - Tj1 && dt < Ta )
        {
            v = v_lim - j_limit * ( Ta - dt ) * ( Ta - dt ) / 2;
        }
        // 2.匀速段
        // a) 匀速度
        if ( dt >= Ta && dt < Ta + Tv )
        {
            v = v_lim;
        }
        // 3.减速段
        // a) 加速度减小
        if ( dt >= T - Td && dt < T - Td + Tj2 )
        {
            v = v_lim + jmin * ( dt - T + Td ) * ( dt - T + Td ) / 2.0;
        }
        // b) 匀加速度
        if ( dt >= T - Td + Tj2 && dt < T - Tj2 )
        {
            v = v_lim + a_limd * ( dt - T + Td - Tj2 / 2 );
        }
        // c) 加速度增大
        if ( dt >= T - Tj2 )
        {
            v = vf - jmin * ( T - dt ) * ( T - dt ) / 2.0;
        }

        if ( dt >= T )
        {
            v = vf;
        }

        v = v * _direction;

        return v;
    }

    double DoubleS::acc( double t ) const
    {
        double a;

        double dt = ( t - _t[ 0 ] );

        double Tj1, Tj2, Ta, Tv, Td, T, v_lim, a_lima, a_limd, j_limit, jmin, p0, pf, v0, vf;

        j_limit = _j[ 1 ];
        //jmin = -j_limit;
        jmin = -_j[ 7 ];

        a_lima = _a[ 1 ];
        a_limd = _a[ 5 ];

        v_lim = _v[ 3 ];

        Tj1 = _t[ 1 ] - _t[ 0 ];
        Ta  = _t[ 3 ] - _t[ 0 ];
        Tv  = _t[ 4 ] - _t[ 3 ];
        T   = _t[ 7 ] - _t[ 0 ];
        Td  = T - Tv - Ta;
        Tj2 = a_limd / jmin;

        p0 = _x[ 0 ];
        pf = _x[ 7 ];
        v0 = _v[ 0 ];
        vf = _v[ 7 ];

        // a) 加速度增大
        if ( dt < Tj1 )
        {
            a = j_limit * dt;
        }

        // b) 匀加速度
        if ( dt >= Tj1 && dt < Ta - Tj1 )
        {
            a = a_lima;
        }
        // c) 加速度减小
        if ( dt >= Ta - Tj1 && dt < Ta )
        {
            a = j_limit * ( Ta - dt );
        }
        // 2.匀速段
        // a) 匀速度
        if ( dt >= Ta && dt < Ta + Tv )
        {
            a = 0;
        }
        // 3.减速段
        // a) 加速度减小
        if ( dt >= T - Td && dt < T - Td + Tj2 )
        {
            a = jmin * ( dt - T + Td );
        }
        // b) 匀加速度
        if ( dt >= T - Td + Tj2 && dt < T - Tj2 )
        {
            a = a_limd;
        }
        // c) 加速度增大
        if ( dt >= T - Tj2 )
        {
            a = jmin * ( T - dt );
        }

        if ( dt >= T )
        {
            a = 0;
        }

        a = a * _direction;

        return a;
    }

    double DoubleS::jerk( double t ) const
    {
        double j;

        double dt = ( t - _t[ 0 ] );

        double Tj1, Tj2, Ta, Tv, Td, T, j_limit, jmin;

        j_limit = _j[ 1 ];
        jmin    = -j_limit;

        Tj1 = _t[ 1 ] - _t[ 0 ];
        Ta  = _t[ 3 ] - _t[ 0 ];
        Tv  = _t[ 4 ] - _t[ 3 ];
        T   = _t[ 7 ] - _t[ 0 ];
        Td  = T - Tv - Ta;
        Tj2 = -_a[ 5 ] / j_limit;

        // a) 加速度增大
        if ( dt < Tj1 )
        {
            j = j_limit;
        }

        // b) 匀加速度
        if ( dt >= Tj1 && dt < Ta - Tj1 )
        {
            j = 0.0;
        }
        // c) 加速度减小
        if ( dt >= Ta - Tj1 && dt < Ta )
        {
            j = jmin;
        }
        // 2.匀速段
        // a) 匀速度
        if ( dt >= Ta && dt < Ta + Tv )
        {
            j = 0.0;
        }
        // 3.减速段
        // a) 加速度减小
        if ( dt >= T - Td && dt < T - Td + Tj2 )
        {
            j = jmin;
        }
        // b) 匀加速度
        if ( dt >= T - Td + Tj2 && dt < T - Tj2 )
        {
            j = 0.0;
        }
        // c) 加速度增大
        if ( dt >= T - Tj2 )
        {
            j = j_limit;
        }

        if ( dt >= T )
        {
            j = 0;
        }

        j = j * _direction;

        return j;
    }

    double DoubleS::max_vel( ) const
    {
        double v = 0;

        for ( int i = 0; i < 7; i++ )
        {
            v = fmax( fabs( v ), fabs( _v[ i ] ) );
        }

        for ( int i = 0; i < 7; i++ )
        {
            if ( v == fabs( _v[ i ] ) )
            {
                v = _v[ i ] * _direction;
            }
        }
        return v;
    }

    double DoubleS::max_acc( ) const
    {
        double a = 0;
        for ( int i = 0; i < 7; i++ )
        {
            a = fmax( fabs( a ), fabs( _a[ i ] ) );
        }
        for ( int i = 0; i < 7; i++ )
        {
            if ( fabs( a ) == fabs( _a[ i ] ) )
            {
                a = _a[ i ] * _direction;
            }
        }
        return a;
    }

#pragma endregion DOUBLES

// 3次5次多项式
#pragma region thrd5th

    void Polynomial::planLinearProfileT( double t0, double tf, double p0, double pf )
    {
        _plannedProfile = false;

        _pType = polynomial_linear_T;

        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }

        // the linear polynomial profile:
        // p(t) = a*(t-t0) + b
        // p(t0) = b = p0
        // p(tf) = pf

        double a, b;

        a = ( pf - p0 ) / ( tf - t0 );

        b = p0;

        _coefficient[ 0 ] = 0.0;  // t^5
        _coefficient[ 1 ] = 0.0;  // t^4
        _coefficient[ 2 ] = 0;    // t^3
        _coefficient[ 3 ] = 0;    // t^2
        _coefficient[ 4 ] = a;    // t^1
        _coefficient[ 5 ] = b;    // t^0

        _plannedProfile = true;
    }

    void Polynomial::plan3rdProfileT( double t0, double tf, double p0,
                                      double pf, double v0, double vf )
    {
        _plannedProfile = false;

        _pType  = polynomial_3rd_T;
        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _v[ 0 ] = v0;
        _v[ 1 ] = vf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;
        _a[ 0 ] = 0.0;
        _a[ 1 ] = 0.0;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }

        // the 3rd polynomial profile:
        // p(t) = a*(t-t0)^3 + b*(t-t0)^2 + c*(t-t0) +d
        // p(t0) = d = p0
        // p(tf) = pf

        // v(t) = 3*a*(t-t0)^2 + 2*b*(t-t0) + c
        // v(t0) = c = v0
        // v(tf)  = vf

        double a, b, c, d, T;

        T = tf - t0;

        d = p0;
        c = v0;
        b = ( 3.0 * pf - vf * T - 2.0 * v0 * T - 3.0 * p0 ) / ( T * T );
        a = ( pf - b * T * T - v0 * T - p0 ) / ( T * T * T );

        _coefficient[ 0 ] = 0.0;  // t^5
        _coefficient[ 1 ] = 0.0;  // t^4
        _coefficient[ 2 ] = a;    // t^3
        _coefficient[ 3 ] = b;    // t^2
        _coefficient[ 4 ] = c;    // t^1
        _coefficient[ 5 ] = d;    // t^0

        _plannedProfile = true;
    }

    void Polynomial::plan5thProfileT( double t0, double tf, double p0,
                                      double pf, double v0, double vf, double a0, double af )
    {
        _pType  = polynomial_5th_T;
        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _v[ 0 ] = v0;
        _v[ 1 ] = vf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;
        _a[ 0 ] = a0;
        _a[ 1 ] = af;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }

        // 5th polynomial profile
        // p(t) = a*(t-t0)^5 + b*(t-t0)^4 + c*(t-t0)^3 + d*(t-t0)^2 + e*(t-t0) + f
        // p(t0) = p0 = f
        // p(tf) = pf

        // v(t) = 5*a*(t-t0)^4 + 4*b*(t-t0)^3 + 3*c*(t-t0)^2 + 2*d*(t-t0) + e
        // v(t0) = v0 = e
        // v(tf) = vf

        // a(t) = 20*a*(t-t0)^3 + 12*(t-t0)^2 + 6*(t-t0) + 2*d
        // a(t0) = a0 = 2*d
        // a(tf) = af

        double a, b, c, d, e, f, T;

        T = tf - t0;

        f = p0;
        e = v0;
        d = a0 * 0.5;
        c = -( 20.0 * f - 20.0 * pf + 12.0 * e * T + 8.0 * vf * T - af * T * T + 6.0 * d * T * T ) / ( 2.0 * T * T * T );
        b = ( 15.0 * f - 15.0 * pf + 8.0 * e * T + 7.0 * vf * T - 1.0 * af * T * T + 3.0 * d * T * T ) / ( T * T * T * T );
        a = -( 12 * f - 12 * pf + 6 * e * T + 6 * vf * T - af * T * T + 2 * d * T * T ) / ( 2 * T * T * T * T * T );

        _coefficient[ 0 ] = a;  // t^5
        _coefficient[ 1 ] = b;  // t^4
        _coefficient[ 2 ] = c;  // t^3
        _coefficient[ 3 ] = d;  // t^2
        _coefficient[ 4 ] = e;  // t^1
        _coefficient[ 5 ] = f;  // t^0

        _plannedProfile = true;
    }

    double Polynomial::scaleToDuration( double newDuration )
    {
        if ( _pType == polynomial_3rd_T )
        {
            plan3rdProfileT( _t[ 0 ], _t[ 0 ] + newDuration, _x[ 0 ], _x[ 1 ], _v[ 0 ], _v[ 1 ] );
        }

        if ( _pType == polynomial_5th_T )
        {
            plan5thProfileT( _t[ 0 ], _t[ 0 ] + newDuration, _x[ 0 ], _x[ 1 ], _v[ 0 ], _v[ 1 ], _a[ 0 ], _a[ 1 ] );
        }

        return getDuration( );
    }

    bool Polynomial::isValidMovement( ) const
    {
        return _plannedProfile;
    }

    double Polynomial::getDuration( ) const
    {
        return _t[ 1 ] - _t[ 0 ];
    }

    double Polynomial::pos( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double p = _coefficient[ 0 ] * t * t * t * t * t + _coefficient[ 1 ] * t * t * t * t + _coefficient[ 2 ] * t * t * t + _coefficient[ 3 ] * t * t + _coefficient[ 4 ] * t + _coefficient[ 5 ];
        return p;
    }

    double Polynomial::vel( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double v = 5 * _coefficient[ 0 ] * t * t * t * t + 4 * _coefficient[ 1 ] * t * t * t + 3 * _coefficient[ 2 ] * t * t + 2 * _coefficient[ 3 ] * t + _coefficient[ 4 ];
        return v;
    }

    double Polynomial::acc( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double a = 20 * _coefficient[ 0 ] * t * t * t + 12 * _coefficient[ 1 ] * t * t + 6 * _coefficient[ 2 ] * t + 2 * _coefficient[ 3 ];
        return a;
    }

    double Polynomial::jerk( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double j = 60 * _coefficient[ 0 ] * t * t + 24 * _coefficient[ 1 ] * t + 6 * _coefficient[ 2 ];
        return j;
    }

    double Polynomial::max_vel( ) const
    {
        double temp, vmax;
        switch ( _pType )
        {
            case polynomial_3rd_T:

                temp = -2.0 * _coefficient[ 3 ] / ( 6.0 * _coefficient[ 2 ] );
                vmax = fmax( fmax( fabs( _v[ 0 ] ), fabs( _v[ 1 ] ) ), fabs( vel( temp ) ) );

                if ( vmax == fabs( vel( _t[ 0 ] ) ) )
                {
                    vmax = vel( _t[ 0 ] );
                }

                if ( vmax == fabs( vel( _t[ 1 ] ) ) )
                {
                    vmax = vel( _t[ 1 ] );
                }

                if ( vmax == fabs( vel( temp ) ) )
                {
                    vmax = vel( temp );
                }
                return vmax;

                break;
            case polynomial_5th_T:

                double a = 20 * _coefficient[ 0 ];
                double b = 12 * _coefficient[ 1 ];
                double c = 6 * _coefficient[ 2 ];
                double d = 2 * _coefficient[ 3 ];

                // get the roots of a*x^3 + b*x^2 + c*x +d = 0;
                double alpha = ( -b * b * b / ( 27.0 * a * a * a ) + -d / ( 2.0 * a ) + b * c / ( 6.0 * a * a ) );
                double beta  = ( c / ( 3.0 * a ) - b * b / ( 9.0 * a * a ) );
                double t1, t2, t3;
                temp = alpha * alpha + beta * beta * beta;
                if ( temp <= 0 )
                {
                    temp = acos( alpha / ( pow( -beta, 1.5 ) ) );
                    t1   = -b / ( 3.0 * a ) + 2.0 * sqrt( -beta ) * cos( temp / 3.0 );
                    t2   = -b / ( 3.0 * a ) + 2.0 * sqrt( -beta ) * cos( ( temp + 2 * PI ) / 3.0 );
                    t3   = -b / ( 3.0 * a ) + 2.0 * sqrt( -beta ) * cos( ( temp - 2 * PI ) / 3.0 );
                }
                else
                {
                    t1 = -b / ( 3.0 * a ) +
                         pow( ( b * c / ( 6.0 * a * a ) - b * b * b / ( 27.0 * a * a * a ) - d / ( 2.0 * a ) + sqrt( temp ) ),
                              1.0 / 3.0 );
                    t1 = t1 + pow( ( b * c / ( 6.0 * a * a ) - b * b * b / ( 27.0 * a * a * a ) - d / ( 2.0 * a ) - sqrt( temp ) ),
                                   1.0 / 3.0 );
                    t2 = t1;
                    t3 = t1;
                }

                vmax = fmax( fmax( fabs( vel( t1 ) ), fabs( vel( t2 ) ) ), fabs( vel( t3 ) ) );
                vmax = fmax( fmax( fabs( _v[ 0 ] ), fabs( _v[ 1 ] ) ), vmax );

                if ( vmax == fabs( vel( t1 ) ) )
                {
                    vmax = vel( t1 );
                }

                if ( vmax == fabs( vel( t2 ) ) )
                {
                    vmax = vel( t2 );
                }

                if ( vmax == fabs( vel( t3 ) ) )
                {
                    vmax = vel( t3 );
                }

                if ( vmax == fabs( vel( _t[ 0 ] ) ) )
                {
                    vmax = vel( _t[ 0 ] );
                }

                if ( vmax == fabs( vel( _t[ 1 ] ) ) )
                {
                    vmax = vel( _t[ 1 ] );
                }

                return vmax;
                break;
        }

        return 0xffffffff;
    }

    double Polynomial::max_acc( ) const
    {
        switch ( _pType )
        {
            case polynomial_3rd_T:

                return fmax( fabs( acc( _t[ 0 ] ) ), fabs( acc( _t[ 1 ] ) ) );

                break;
            case polynomial_5th_T:

                double a = 60 * _coefficient[ 0 ];
                double b = 24 * _coefficient[ 1 ];
                double c = 6 * _coefficient[ 2 ];

                // get the roots of a*x^2 + b*x + c = 0;
                double temp = b * b - 4 * a * c, amax;
                if ( temp >= 0 )
                {
                    double t1 = ( -b + sqrt( temp ) ) / ( 2 * a );
                    double t2 = ( -b - sqrt( temp ) ) / ( 2 * a );
                    amax      = fmax( fabs( acc( _t[ 0 ] ) ), fabs( acc( _t[ 1 ] ) ) );
                    amax      = fmax( amax, fmax( fabs( acc( t1 ) ), fabs( acc( t2 ) ) ) );
                    if ( amax == fabs( acc( t1 ) ) )
                    {
                        amax = acc( t1 );
                    }
                    if ( amax == fabs( acc( t2 ) ) )
                    {
                        amax = acc( t2 );
                    }
                }
                else
                {
                    amax = fmax( fabs( acc( _t[ 0 ] ) ), fabs( acc( _t[ 1 ] ) ) );
                }
                if ( amax == fabs( acc( _t[ 0 ] ) ) )
                {
                    amax = acc( _t[ 0 ] );
                }
                if ( amax == fabs( acc( _t[ 1 ] ) ) )
                {
                    amax = acc( _t[ 1 ] );
                }

                return amax;
                break;
        }

        return 0;
    }

    double Polynomial::max_jerk( ) const
    {
        switch ( _pType )
        {
            case polynomial_3rd_T:
                return jerk( _t[ 0 ] );
                break;
            case polynomial_5th_T:
                double a = 120 * _coefficient[ 0 ];
                double b = 24 * _coefficient[ 1 ];
                double t = -b / ( a * 2.0 );

                double jmax = fmax( fabs( jerk( _t[ 0 ] ) ), fabs( jerk( _t[ 1 ] ) ) );
                jmax        = fmax( jmax, fabs( jerk( t ) ) );
                if ( jmax == fabs( jerk( t ) ) )
                {
                    jmax = jerk( t );
                }

                if ( jmax == fabs( jerk( _t[ 0 ] ) ) )
                {
                    jmax = jerk( _t[ 0 ] );
                }
                if ( jmax == fabs( jerk( _t[ 1 ] ) ) )
                {
                    jmax = jerk( _t[ 1 ] );
                }

                return jmax;

                break;
        }

        return 0;
    }

#pragma endregion thrd5th

// 三角谐波
#pragma region HARMONIC

    void Trigonometric::planHarmonicProfileT( double p0,
                                              double pf, double t0, double tf )
    {
        _pType  = harmonic_T;
        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }
        _plannedProfile = true;
    }

    void Trigonometric::planCycloidalProfileT( double p0,
                                               double pf, double t0, double tf )
    {
        _pType  = cycloidal_T;
        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }
        _plannedProfile = true;
    }

    void Trigonometric::planEllipticProfileT( double p0,
                                              double pf, double t0, double tf, double n )
    {
        _pType  = harmonic_T;
        _t[ 0 ] = t0;
        _t[ 1 ] = tf;
        _x[ 0 ] = p0;
        _x[ 1 ] = pf;
        _n      = n;

        if ( tf <= t0 )
        {
            _plannedProfile = false;
            //std::cerr<<("start and terminal time cannot be the same!\n");
            return;
        }
        _plannedProfile = true;
    }

    double Trigonometric::scaleToDuration( double newDuration )
    {
        _t[ 1 ] = _t[ 0 ] + newDuration;
        return getDuration( );
    }

    bool Trigonometric::isValidMovement( ) const
    {
        return _plannedProfile;
    }

    double Trigonometric::getDuration( ) const
    {
        return _t[ 1 ] - _t[ 0 ];
    }

    double Trigonometric::pos( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double p = 0.0;

        switch ( _pType )
        {
            case harmonic_T:
                p = ( 0.5 * ( 1 - cos( M_PI * t ) ) );
                p = _x[ 0 ] + ( _x[ 1 ] - _x[ 0 ] ) * p;
                break;
            case cycloidal_T:
                p = t - 1.0 / ( 2 * M_PI ) * sin( 2 * M_PI * t );
                p = _x[ 0 ] + ( _x[ 1 ] - _x[ 0 ] ) * p;
                break;
            case elliptic_T:
                double alpha = ( _n * _n - 1.0 ) / ( _n * _n );
                p            = 0.5 * ( 1 - cos( M_PI * t ) / sqrt( 1.0 - alpha * ( sin( M_PI * t ) * sin( M_PI * t ) ) ) );
                p            = _x[ 0 ] + ( _x[ 1 ] - _x[ 0 ] ) * p;
                break;
        }

        return p;
    }

    double Trigonometric::vel( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double v = 0.0;

        switch ( _pType )
        {
            case harmonic_T:
                v = M_PI / 2.0 * sin( M_PI * t );
                v = v * ( _x[ 1 ] - _x[ 0 ] ) / ( _t[ 1 ] - _t[ 0 ] );
                break;
            case cycloidal_T:
                v = 1.0 - cos( 2.0 * M_PI * t );
                v = v * ( _x[ 1 ] - _x[ 0 ] ) / ( _t[ 1 ] - _t[ 0 ] );
                break;
            case elliptic_T:
                double alpha = ( _n * _n - 1.0 ) / ( _n * _n );
                double temp  = sin( M_PI * t ) * sin( M_PI * t );
                temp         = 1.0 - alpha * temp;
                temp         = temp * temp * temp;
                v            = M_PI / 2.0 * sin( M_PI * t ) / ( _n * _n ) / sqrt( temp );
                v            = v * ( _x[ 1 ] - _x[ 0 ] ) / ( _t[ 1 ] - _t[ 0 ] );
                break;
        }

        return v;
    }

    double Trigonometric::acc( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double a = 0.0;
        double T = ( ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) );

        switch ( _pType )
        {
            case harmonic_T:
                a = M_PI * M_PI / 2.0 * cos( M_PI * t );
                a = a * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case cycloidal_T:
                a = 2.0 * M_PI * sin( 2.0 * M_PI * t );
                a = a * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case elliptic_T:
                double alpha = ( _n * _n - 1.0 ) / ( _n * _n );
                double x     = sin( M_PI * t ) * sin( M_PI * t );
                double y     = 1.0 + 2.0 * alpha * x;
                double z     = 1.0 - alpha * x;
                z            = z * z * z * z * z;
                a            = M_PI * M_PI / 2.0 * cos( M_PI * T ) * y / ( _n * _n ) / sqrt( z );
                a            = a * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
        }

        return a;
    }

    double Trigonometric::jerk( double t ) const
    {
        t = ( t - _t[ 0 ] );

        if ( t < 0 )
        {
            t = 0;
        }

        if ( t > ( _t[ 1 ] - _t[ 0 ] ) )
        {
            t = _t[ 1 ] - _t[ 0 ];
        }

        double j = 0.0;
        double T = ( ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) );

        switch ( _pType )
        {
            case harmonic_T:
                j = -( M_PI * M_PI * M_PI ) / 2.0 * sin( M_PI * t );
                j = j * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case cycloidal_T:
                j = 4.0 * M_PI * M_PI * cos( 2 * M_PI * t );
                j = j * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case elliptic_T:
                j = 0;
                break;
        }

        return j;
    }

    double Trigonometric::max_vel( ) const
    {
        double v = 0.0;
        double T = ( ( _t[ 1 ] - _t[ 0 ] ) );

        switch ( _pType )
        {
            case harmonic_T:
                v = PI / 2.0;
                v = v * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case cycloidal_T:
                v = 2.0 * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case elliptic_T:
                v = 2.0 * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
        }

        return v;
    }

    double Trigonometric::max_acc( ) const
    {
        double a = 0.0;
        double T = ( ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) );

        switch ( _pType )
        {
            case harmonic_T:
                a = 0.5 * M_PI * M_PI;
                a = a * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case cycloidal_T:
                a = 2.0 * M_PI;
                a = a * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case elliptic_T:
                a = 4.0 * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
        }

        return a;
    }

    double Trigonometric::max_jerk( ) const
    {
        double j = 0.0;
        double T = ( ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) * ( _t[ 1 ] - _t[ 0 ] ) );

        switch ( _pType )
        {
            case harmonic_T:
                j = 0.5 * M_PI * M_PI * M_PI;
                j = j * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case cycloidal_T:
                j = 4.0 * M_PI * M_PI;
                j = j * ( _x[ 1 ] - _x[ 0 ] ) / T;
                break;
            case elliptic_T:
                j = 0;
                break;
        }

        return j;
    }

#pragma endregion HARMONIC

}  // namespace rocos