#include <Eigen/Dense>
#include <atomic>
#include <chrono>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/frames.hpp>
#include <kdl/jntarray.hpp>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Log.h>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <std_msgs/Int64.h>
#include <std_msgs/Int64MultiArray.h>
#include <thread>
#include <trac_ik/trac_ik.hpp>
#include <ur_rtde/rtde_control_interface.h>
#include <ur_rtde/rtde_io_interface.h>
#include <ur_rtde/rtde_receive_interface.h>
#include <velocity_profile/interpolate.h>

using namespace KDL;
using namespace std;
using namespace ur_rtde;
using namespace std::chrono;

/**
 * @brief 四元素插值公式
 *
 * @param end 终止四元素
 * @param s 百分比
 * @return std::vector< double > 在s位置处的四元素
 */
std::vector< double > UnitQuaternion_intep( const std::vector< double >& start, const std::vector< double >& end, double s )
{
    if ( s > 1 || s < 0 )
    {
        std::cerr << "values of S outside interval [0,1]" << std::endl;
    }

    double cosTheta               = start[ 0 ] * end[ 0 ] + start[ 1 ] * end[ 1 ] + start[ 2 ] * end[ 2 ] + start[ 3 ] * end[ 3 ];
    std::vector< double > start_2 = start;

    //** 这里是为了取最短路径 **//
    if ( cosTheta < 0 )
    {
        for ( int i = 0; i < 4; i++ ) start_2[ i ] *= -1;
        cosTheta *= -1;
    }
    //**-------------------------------**//

    double theta = acos( cosTheta );
    if ( theta == 0 || s == 0 )
        return start_2;
    else
    {
        double coefficient_1 = sin( ( 1 - s ) * theta ) / sin( theta );
        double coefficient_2 = sin( ( s )*theta ) / sin( theta );

        return std::vector< double >{
            coefficient_1 * start_2[ 0 ] + coefficient_2 * end[ 0 ],
            coefficient_1 * start_2[ 1 ] + coefficient_2 * end[ 1 ],
            coefficient_1 * start_2[ 2 ] + coefficient_2 * end[ 2 ],
            coefficient_1 * start_2[ 3 ] + coefficient_2 * end[ 3 ] };
    }
}

/**
 * @brief 直线规划（使用doubleS速度曲线）
 *
 * @param start 起始位姿
 * @param end 终止位姿
 * @return std::vector< KDL::Frame > 轨迹
 */
std::vector< KDL::Frame > moveL( KDL::Frame start, KDL::Frame end, double path_lenth = 0 )
{
    //** 变量初始化 **//
    KDL::Vector Pstart = start.p;
    KDL::Vector Pend   = end.p;
    std::vector< KDL::Frame > traj_1;
    double s = 0;
    rocos::R_INTERP doubleS;
    bool isplanned       = doubleS.planDoubleSPorfile( 0, 0, 1, 0, 0, 0.1 / path_lenth, 1 / path_lenth, 6 / path_lenth );
    const double timegap = 0.002;  // 500hz
    assert( isplanned && ( doubleS.getDuration( ) > 0 ) && "doubleS planning is failed" );
    int N = doubleS.getDuration( ) / timegap;
    std::vector< double > Quaternion_start{ 0, 0, 0, 0 };
    std::vector< double > Quaternion_end{ 0, 0, 0, 0 };
    std::vector< double > Quaternion_interp{ 0, 0, 0, 0 };
    start.M.GetQuaternion( Quaternion_start.at( 0 ), Quaternion_start.at( 1 ), Quaternion_start.at( 2 ), Quaternion_start.at( 3 ) );
    end.M.GetQuaternion( Quaternion_end.at( 0 ), Quaternion_end.at( 1 ), Quaternion_end.at( 2 ), Quaternion_end.at( 3 ) );

    //**-------------------------------**//

    //** 轨迹计算 **//
    for ( int i = 0; i <= N; i++ )
    {
        s                 = doubleS.pos( timegap * i );
        KDL::Vector P     = Pstart + ( Pend - Pstart ) * s;
        Quaternion_interp = UnitQuaternion_intep( Quaternion_start, Quaternion_end, s );
        traj_1.push_back( KDL::Frame( KDL::Rotation::Quaternion( Quaternion_interp[ 0 ], Quaternion_interp[ 1 ], Quaternion_interp[ 2 ], Quaternion_interp[ 3 ] ), P ) );
    }
    //**-------------------------------**//
    return traj_1;
}


/**
 * @brief 计算直线轨迹的速度向量
 * 
 * @param Cartesian_vel 计算的笛卡尔速度向量
 * @param t 求解时间点
 * @param start 起始位姿
 * @param end 终止位姿
 * @return int 
 */
int moveL_vel( KDL::Twist& Cartesian_vel, double t, KDL::Frame start, KDL::Frame end )
{
    KDL::Vector vel   = end.p - start.p;
    double path_lenth = vel.Normalize( );

    rocos::R_INTERP doubleS;
    bool isplanned    = doubleS.planDoubleSPorfile( 0, 0, 1, 0, 0, 0.1 / path_lenth, 1 / path_lenth, 6 / path_lenth );
    vel               = vel * doubleS.vel( t );
    Cartesian_vel.vel = vel;

    KDL::Rotation rot_start_end = start.M.Inverse( ) * end.M;
    KDL::Vector aixs            = rot_start_end.GetRot( );
    aixs.Normalize( );
    Cartesian_vel.rot = aixs * doubleS.vel( t );

    return 0;
}

class ur_force_control
{
public:
    static const unsigned int joint_num{ 6 };
    //机器人读取变量
private:
    KDL::JntArray max_joints_vel{ joint_num };

    std::vector< double > TCP_force{ 0, 0, 0 };
    std::vector< double > force_pos_offset{ 0, 0, 0 };
    std::vector< double > force_vel_offset{ 0, 0, 0 };
    std::vector< double > force_last_vel_offset{ 0, 0, 0 };
    std::vector< double > force_acc_offset{ 0, 0, 0 };
    std::vector< double > force_last_acc_offset{ 0, 0, 0 };

    std::vector< double > TCP_torque{ 0, 0, 0 };
    std::vector< double > torque_pos_offset{ 0, 0, 0 };
    std::vector< double > torque_vel_offset{ 0, 0, 0 };
    std::vector< double > torque_last_vel_offset{ 0, 0, 0 };
    std::vector< double > torque_acc_offset{ 0, 0, 0 };
    std::vector< double > torque_last_acc_offset{ 0, 0, 0 };

    std::vector< double > M{ 30, 30, 30, 30, 30, 30 };
    std::vector< double > K{ 150, 150, 150, 100, 100, 100 };
    std::vector< double > B{ 30, 30, 30, 30, 30, 30 };

    KDL::Twist _Cartesian_vel;
    double dt{ 0.002 };

    std::shared_ptr< RTDEControlInterface > _ur_control_ptr{ nullptr };

public:
    ur_force_control( std::shared_ptr< RTDEControlInterface > ur_control_ptr )
        : _ur_control_ptr( ur_control_ptr )
    {
        for ( int i{ 0 }; i < joint_num; i++ )
            B[ i ] = 2 * 1 * sqrt( M[ i ] * K[ i ] );

        for ( int i{ 0 }; i < joint_num; i++ )
            max_joints_vel( i ) = 40 * M_PI / 180;
    }

    /**
     * @brief
     *
     * @param pos_offset  当前位置-期望位置 -  轨迹跟踪引起的误差 = 力引起的误差
     * @param vel_offset 当前速度 -  期望速度 -  轨迹跟踪引起的误差=  力引起的误差
     */
    void init( const KDL::Vector& pos_offset, const KDL::Vector& vel_offset )
    {
        for ( int i{ 0 }; i < 3; i++ )
        {
            force_pos_offset[ i ]  = pos_offset[ i ];
            force_vel_offset[ i ]  = vel_offset[ i ];
            torque_pos_offset[ i ] = pos_offset[ 3 + i ];
            torque_vel_offset[ i ] = vel_offset[ 3 + i ];
        }
    }

    void calculate_translate( )
    {
        for ( int i{ 0 }; i < 3; i++ )
        {
            force_acc_offset[ i ] = ( TCP_force[ i ] - B[ i ] * force_vel_offset[ i ] - K[ i ] * force_pos_offset[ i ] ) / M[ i ];
            force_vel_offset[ i ] = dt * ( force_acc_offset[ i ] + force_last_acc_offset[ i ] ) / 2 + force_vel_offset[ i ];
            force_pos_offset[ i ] = dt * ( force_vel_offset[ i ] + force_last_vel_offset[ i ] ) / 2 + force_pos_offset[ i ];

            force_last_acc_offset[ i ] = force_acc_offset[ i ];
            force_last_vel_offset[ i ] = force_vel_offset[ i ];

            _Cartesian_vel.vel[ i ] = force_vel_offset[ i ];
        }
    }

    KDL::Rotation calculate_rotation( )
    {
        KDL::Vector delta_rot;
        KDL::Vector current_rot;
        KDL::Rotation template_rot;

        for ( int i{ 0 }; i < 3; i++ )
        {
            torque_acc_offset[ i ] = ( TCP_torque[ i ] - B[ i ] * torque_vel_offset[ i ] - K[ i ] * torque_pos_offset[ i ] ) / M[ i ];
            torque_vel_offset[ i ] = dt * ( torque_acc_offset[ i ] + torque_last_acc_offset[ i ] ) / 2 + torque_vel_offset[ i ];

            delta_rot( i )   = dt * ( torque_vel_offset[ i ] + torque_last_vel_offset[ i ] ) / 2;
            current_rot( i ) = torque_pos_offset[ i ];

            torque_last_acc_offset[ i ] = torque_acc_offset[ i ];
            torque_last_vel_offset[ i ] = torque_vel_offset[ i ];

            _Cartesian_vel.rot[ i ] = torque_vel_offset[ i ];
        }

        template_rot = KDL::Rotation::Rot( delta_rot, delta_rot.Norm( ) ) * KDL::Rotation::Rot( current_rot, current_rot.Norm( ) );

        current_rot = template_rot.GetRot( );

        for ( int i{ 0 }; i < 3; i++ )
            torque_pos_offset[ i ] = current_rot[ i ];

        return template_rot;
    }

    int calculate( KDL::Frame& pos_offset, KDL::Twist& Cartesian_vel )
    {
        calculate_translate( );

        for ( int i{ 0 }; i < 3; i++ )
        {
            pos_offset.p[ i ] = force_pos_offset[ i ];
        }

        pos_offset.M = calculate_rotation( );

        Cartesian_vel = _Cartesian_vel;

        return 0;
    }

    void check_before_move( KDL::JntArray& joints_vel )
    {
        for ( int i{ 0 }; i < joint_num; i++ )
        {
            if ( joints_vel( i ) > max_joints_vel( i ) || TCP_force[ i > 2 ? 2 : i ] > 70 || TCP_torque[ i > 2 ? 2 : i ] > 30 )
            {
                PLOG_ERROR << "joint [" << i << "]  速度过快";
                PLOG_ERROR << "目标速度  = " << joints_vel( i );
                PLOG_ERROR << "允许最大速度  = " << max_joints_vel( i );

                if ( _ur_control_ptr )
                {
                    _ur_control_ptr->servoStop( );
                    _ur_control_ptr->stopScript( );
                }
                exit( 0 );
            }
        }
    }

    void set_force( double force_x, double force_y, double force_z )
    {
        TCP_force[ 0 ] = force_x;
        TCP_force[ 1 ] = force_y;
        TCP_force[ 2 ] = force_z;
        // PLOG_DEBUG.printf( "TCP_force  = %f %f %f", TCP_force[ 0 ], TCP_force[ 1 ], TCP_force[ 2 ] );
    }

    void set_torque( double tor_que_x, double tor_que_y, double tor_que_z )
    {
        TCP_torque[ 0 ] = tor_que_x;
        TCP_torque[ 1 ] = tor_que_y;
        TCP_torque[ 2 ] = tor_que_z;
        // PLOG_DEBUG.printf( "TCP_torque  = %f %f %f", TCP_torque[ 0 ], TCP_torque[ 1 ], TCP_torque[ 2 ] );
    }

    void set_damp( double value )
    {
        static double damp = 1;
        damp += value;

        for ( int i{ 0 }; i < joint_num; i++ )
            B[ i ] = 2 * damp * sqrt( M[ i ] * K[ i ] );

        PLOG_DEBUG << "damp  = " << damp;
    }
};

void RunTeleop_damp( const std_msgs::Int64::ConstPtr& msg, ur_force_control* ur_admittance_ptr )
{
    if ( msg->data == 9 )
        ur_admittance_ptr->set_damp( 0.1 );
    else if ( msg->data == 10 )
        ur_admittance_ptr->set_damp( -0.1 );
}

void RunTeleop_force_torque( const std_msgs::Int64MultiArray::ConstPtr& msg, ur_force_control* ur_admittance_ptr )
{
    ur_admittance_ptr->set_force( msg->data[ 0 ], msg->data[ 1 ], msg->data[ 2 ] );
    ur_admittance_ptr->set_torque( msg->data[ 3 ], msg->data[ 4 ], msg->data[ 5 ] );
}

int main( int argc, char* argv[] )
{
    //** 变量初始化 **//
    static plog::ColorConsoleAppender< plog::TxtFormatter > consoleAppender;

    plog::init( plog::debug, &consoleAppender );  // Initialize the logger.
    ros::init( argc, argv, "moveL_admittance_Control" );
    ros::NodeHandle nh;
    KDL::Chain UR_chain;
    TRAC_IK::TRAC_IK tracik_solver( "base", "wrist_3_link" );
    std::shared_ptr< RTDEControlInterface > rtde_control_ptr{ nullptr };
    std::shared_ptr< RTDEReceiveInterface > rtde_receive_ptr{ nullptr };

    tracik_solver.getKDLChain( UR_chain );
    ChainFkSolverPos_recursive ur_fk{ UR_chain };
    KDL::ChainIkSolverVel_pinv ik_vel{ UR_chain };

    bool is_remoted_control{ false };
    nh.param( "/moveL_admittance_Control/is_remoted_control", is_remoted_control, false );
    if ( is_remoted_control )
    {
        rtde_control_ptr.reset( new RTDEControlInterface{ "192.168.3.101", RTDEControlInterface::FLAG_USE_EXT_UR_CAP } );
        rtde_receive_ptr.reset( new RTDEReceiveInterface{ "192.168.3.101" } );
    }

    ur_force_control ur_admittance{ rtde_control_ptr };
    ros::Subscriber sub        = nh.subscribe< std_msgs::Int64 >( "dragging_teleop", 1000, boost::bind( &RunTeleop_damp, _1, &ur_admittance ) );
    ros::Publisher publisher   = nh.advertise< sensor_msgs ::JointState >( "my_joint_state", 1000 );
    ros::Subscriber sub_torque = nh.subscribe< std_msgs::Int64MultiArray >( "force_torque_control", 1000, boost::bind( &RunTeleop_force_torque, _1, &ur_admittance ) );

    ros::AsyncSpinner spinner( 2 );
    spinner.start( );

    const std::vector< double > test_joint{ 0 * M_PI / 180, ( -145 + 30 ) * M_PI / 180, ( -80 - 30 ) * M_PI / 180, -45 * M_PI / 180, 90 * M_PI / 180, 0 };
    KDL::JntArray joint_init{ 6 };
    KDL::JntArray joint_result{ 6 };

    sensor_msgs::JointState command_msg;
    KDL::JntArray real_joint( 6 );

    //**-------------------------------**//

    command_msg.name = ros::V_string{ "shoulder_pan_joint", "shoulder_lift_joint", "elbow_joint", "wrist_1_joint", "wrist_2_joint", "wrist_3_joint" };
    command_msg.position.resize( ur_force_control::joint_num );

    if ( is_remoted_control )
    {
        rtde_control_ptr->moveJ( test_joint, 0.3, 0.2, false );  //起始位置初始化}
    }

    for ( int i = 0; i < ur_force_control::joint_num; i++ )
    {
        joint_init.data[ i ]   = test_joint[ i ];
        joint_result.data[ i ] = test_joint[ i ];
        real_joint( i )        = test_joint[ i ];
    }

    KDL::Frame frame_init;
    KDL::Frame frame_target;
    KDL::Frame frame_offset;
    KDL::JntArray joints_vel( 6 );
    KDL::Twist admittance_vel;
    KDL::Twist traj_vel;
    KDL::Twist Cartesian_vel;
    std::vector< double > tcp_zero_force_torque( 6, 0 );
    std::vector< double > tcp_force_torque( 6, 0 );
    std::vector< double > joint_command( 6, 0 );

    if ( is_remoted_control )
    {
        sleep( 3 );
        tcp_zero_force_torque = rtde_receive_ptr->getActualTCPForce( );
        PLOG_DEBUG << tcp_zero_force_torque[ 0 ];
        PLOG_DEBUG << tcp_zero_force_torque[ 1 ];
        PLOG_DEBUG << tcp_zero_force_torque[ 2 ];
    }

    //** 轨迹计算**//
    KDL::Frame start;
    KDL::Vector path_lenth = KDL::Vector( 0, -0.3, 0 );
    ur_fk.JntToCart( joint_init, start );
    KDL::Frame end                   = start * KDL::Frame( path_lenth );
    std::vector< KDL::Frame > traj_1 = moveL( start, end, path_lenth.Norm( ) );
    std::vector< KDL::Frame > traj_2 = moveL( end, start, path_lenth.Norm( ) );

    // //**-------------------------------**//

    int count{ 0 };

    while ( ros::ok( ) )
    {
        auto t_start = steady_clock::now( );

        ur_admittance.calculate( frame_offset, admittance_vel );

        if ( count >= ( traj_1.size( ) + traj_2.size( ) ) )
        {
            count = 0;
        }

        if ( count < traj_1.size( ) )//直线向前
            frame_target = frame_offset * traj_1[ count ];
        else if ( count < ( traj_1.size( ) + traj_2.size( ) ) )//直线向后
            frame_target = frame_offset * traj_2[ count - traj_1.size( ) ];

        int success = tracik_solver.CartToJnt( joint_init, frame_target, joint_result );
        //!IK失败就用笛卡尔速度转关节速度，再计算关节命令
        if ( success < 0 )
        {
            PLOG_ERROR << "****CartToJnt fail  ";

            if ( count < traj_1.size( ) )
                moveL_vel( traj_vel, count * 0.002, start, end );
            else if ( count < ( traj_1.size( ) + traj_2.size( ) ) )
                moveL_vel( traj_vel, ( count - traj_1.size( ) ) * 0.002, end, start );

            Cartesian_vel.vel = traj_vel.vel + admittance_vel.vel;
            Cartesian_vel.rot = traj_vel.rot + admittance_vel.rot;
            ik_vel.CartToJnt( joint_init, Cartesian_vel, joints_vel );
            KDL::Multiply( joints_vel, 0.002, joints_vel );
            KDL::Add( joint_init, joints_vel, joint_result );
        }
        //IK成功
        else
            joint_init = joint_result;

        for ( int i = 0; i < ur_force_control::joint_num; i++ )
        {
            command_msg.position[ i ] = joint_result( i );
        }

        publisher.publish( command_msg );

        if ( is_remoted_control )
        {
            double velocity       = 0.5;
            double acceleration   = 0.5;
            double dt             = 1.0 / 500;  // 2ms
            double lookahead_time = 0.1;
            double gain           = 300;

            tcp_force_torque = rtde_receive_ptr->getActualTCPForce( );

            for ( int i = 0; i < ur_force_control::joint_num; i++ )
            {
                joint_command[ i ] = joint_result( i );
                tcp_force_torque[ i ] -= tcp_zero_force_torque[ i ];
            }

            rtde_control_ptr->servoJ( joint_command, velocity, acceleration, dt, lookahead_time, gain );

            ur_admittance.set_force( tcp_force_torque[ 0 ], tcp_force_torque[ 1 ], tcp_force_torque[ 2 ] );
            ur_admittance.set_torque( tcp_force_torque[ 3 ], tcp_force_torque[ 4 ], tcp_force_torque[ 5 ] );
        }

        auto t_stop     = steady_clock::now( );
        auto t_duration = std::chrono::duration< double >( t_stop - t_start );
        if ( t_duration.count( ) < 0.002 )
        {
            std::this_thread::sleep_for( std::chrono::duration< double >( 0.002 - t_duration.count( ) ) );
        }

        count++;
    }

    if ( is_remoted_control )
    {
        rtde_control_ptr->servoStop( );
        rtde_control_ptr->stopScript( );
    }

    return 0;
}
