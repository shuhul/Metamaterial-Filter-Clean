// Privileged, deliberately small protocol: Status, List, Stop.
// Installed executable/configuration are administrator-owned, never workspace-loaded.
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.IO.Pipes;
using System.Runtime.InteropServices;
using System.Security.AccessControl;
using System.Security.Cryptography;
using System.Security.Principal;
using System.Text;
using System.Web.Script.Serialization;

class SonnetEmBroker
{
    const uint Query = 0x1000, Terminate = 1, Synchronize = 0x100000;
    [DllImport("kernel32.dll", SetLastError=true)] static extern IntPtr OpenProcess(uint access, bool inherit, int id);
    [DllImport("kernel32.dll")] static extern bool CloseHandle(IntPtr handle);
    [DllImport("kernel32.dll", SetLastError=true)] static extern bool QueryFullProcessImageName(IntPtr handle, uint flags, StringBuilder path, ref int length);
    [DllImport("kernel32.dll", SetLastError=true)] static extern bool GetProcessTimes(IntPtr handle, out long creation, out long exit, out long kernel, out long user);
    [DllImport("kernel32.dll", SetLastError=true)] static extern bool ProcessIdToSessionId(int id, out uint session);
    [DllImport("kernel32.dll", SetLastError=true)] static extern bool TerminateProcess(IntPtr handle, uint code);
    [DllImport("kernel32.dll")] static extern uint WaitForSingleObject(IntPtr handle, uint milliseconds);
    static JavaScriptSerializer Json = new JavaScriptSerializer { MaxJsonLength=16384, RecursionLimit=12 };
    static string Folder = AppDomain.CurrentDomain.BaseDirectory;
    static string AllowedImage, AllowedHash, PipeName, UserSid;

    static string Text(Dictionary<string,object> values, string key)
    { return values.ContainsKey(key) ? Convert.ToString(values[key], CultureInfo.InvariantCulture) : ""; }
    static Dictionary<string,object> Reply(string status)
    { return new Dictionary<string,object> { {"Status", status}, {"Utc",DateTime.UtcNow.ToString("o")} }; }
    static void Audit(object data)
    {
        string path=Path.Combine(Folder,"audit.jsonl");
        if(File.Exists(path) && new FileInfo(path).Length>5*1024*1024)
        {
            string previous=path+".previous"; if(File.Exists(previous)) File.Delete(previous);
            File.Move(path,previous);
        }
        File.AppendAllText(path,Json.Serialize(data)+Environment.NewLine,Encoding.UTF8);
    }
    static string Image(IntPtr handle)
    {
        StringBuilder buffer=new StringBuilder(32768); int length=buffer.Capacity;
        if (!QueryFullProcessImageName(handle,0,buffer,ref length)) throw new System.ComponentModel.Win32Exception();
        return Path.GetFullPath(buffer.ToString());
    }
    static long Creation(IntPtr handle)
    {
        long creation,exit,kernel,user;
        if (!GetProcessTimes(handle,out creation,out exit,out kernel,out user)) throw new System.ComponentModel.Win32Exception();
        return creation;
    }
    static bool IsAllowed(string path)
    { return String.Equals(path,AllowedImage,StringComparison.OrdinalIgnoreCase); }
    static bool HashMatches()
    {
        using (SHA256 hash=SHA256.Create()) using (FileStream file=File.OpenRead(AllowedImage))
            return String.Equals(BitConverter.ToString(hash.ComputeHash(file)).Replace("-",""),AllowedHash,StringComparison.OrdinalIgnoreCase);
    }
    static object List()
    {
        List<object> found=new List<object>();
        foreach (Process process in Process.GetProcessesByName("em")) using(process)
        {
            IntPtr handle=OpenProcess(Query,false,process.Id);
            if (handle==IntPtr.Zero) continue;
            try
            {
                string path=Image(handle); if (!IsAllowed(path)) continue;
                uint session; if (!ProcessIdToSessionId(process.Id,out session)) continue;
                found.Add(new Dictionary<string,object> { {"ProcessId",process.Id}, {"Name","em.exe"},
                    {"SessionId",session}, {"CreationUtc",DateTime.FromFileTimeUtc(Creation(handle)).ToString("o")},
                    {"ExecutablePath",path} });
            }
            catch (System.ComponentModel.Win32Exception) { }
            finally { CloseHandle(handle); }
        }
        Dictionary<string,object> result=Reply("OK"); result["Processes"]=found; return result;
    }
    static object Stop(Dictionary<string,object> request)
    {
        object raw;
        if (!request.TryGetValue("Target",out raw)) return Reply("InvalidTarget");
        Dictionary<string,object> target=raw as Dictionary<string,object>;
        if (target==null || Text(target,"Name")!="em.exe") return Reply("InvalidTarget");
        int id; uint expectedSession; DateTime expected;
        if (!Int32.TryParse(Text(target,"ProcessId"),out id) || id<=0 ||
            !UInt32.TryParse(Text(target,"SessionId"),out expectedSession) ||
            !DateTime.TryParse(Text(target,"CreationUtc"),CultureInfo.InvariantCulture,
                DateTimeStyles.RoundtripKind,out expected)) return Reply("InvalidTarget");
        Dictionary<string,object> result=Reply("Failed"); result["ProcessId"]=id;
        // Retaining this handle prevents PID reuse between validation and termination.
        IntPtr handle=OpenProcess(Query|Terminate|Synchronize,false,id);
        if (handle==IntPtr.Zero)
        {
            int error=Marshal.GetLastWin32Error();
            result["Status"]=error==87 ? "AlreadyExited" : "AccessDenied"; result["WindowsError"]=error;
            return result;
        }
        try
        {
            string path=Image(handle); uint session;
            if (!IsAllowed(path)) { result["Status"]="NotAllowedExecutable"; return result; }
            if (!ProcessIdToSessionId(id,out session) || session!=expectedSession ||
                Math.Abs((DateTime.FromFileTimeUtc(Creation(handle))-expected.ToUniversalTime()).TotalMilliseconds)>=1)
            { result["Status"]="IdentityChanged"; return result; }
            if (!HashMatches()) { result["Status"]="InstalledBinaryChanged"; return result; }
            if (!TerminateProcess(handle,1)) throw new System.ComponentModel.Win32Exception();
            result["Status"]=WaitForSingleObject(handle,5000)==0 ? "Stopped" : "TerminationRequested";
            return result;
        }
        finally { CloseHandle(handle); }
    }
    static object Handle(Dictionary<string,object> request)
    {
        string action=Text(request,"Action");
        if (action=="Status")
        {
            Dictionary<string,object> status=Reply("Ready");
            status["Identity"]=WindowsIdentity.GetCurrent().Name;
            status["BrokerProcessId"]=Process.GetCurrentProcess().Id;
            status["AllowedImage"]=AllowedImage; status["AllowedUserSid"]=UserSid;
            status["Version"]=1; return status;
        }
        if (action=="List") return List();
        if (action=="Stop") return Stop(request);
        return Reply("InvalidAction");
    }
    static byte[] ReadBytes(Stream pipe,int count)
    {
        byte[] bytes=new byte[count]; int offset=0;
        while(offset<count)
        {
            IAsyncResult pending=pipe.BeginRead(bytes,offset,count-offset,null,null);
            if(!pending.AsyncWaitHandle.WaitOne(3000)) throw new IOException("Request read timed out");
            int n=pipe.EndRead(pending); pending.AsyncWaitHandle.Close();
            if(n==0) throw new EndOfStreamException(); offset+=n;
        }
        return bytes;
    }
    static void WriteBytes(Stream pipe,byte[] bytes)
    {
        IAsyncResult pending=pipe.BeginWrite(bytes,0,bytes.Length,null,null);
        if(!pending.AsyncWaitHandle.WaitOne(3000)) throw new IOException("Response write timed out");
        pipe.EndWrite(pending); pending.AsyncWaitHandle.Close();
    }
    static int Main()
    {
        try
        {
            Dictionary<string,object> config=Json.Deserialize<Dictionary<string,object>>(File.ReadAllText(Path.Combine(Folder,"config.json")));
            AllowedImage=Path.GetFullPath(Text(config,"AllowedImage")); AllowedHash=Text(config,"AllowedHash");
            PipeName=Text(config,"PipeName"); UserSid=Text(config,"UserSid");
            if (WindowsIdentity.GetCurrent().User.Value!="S-1-5-18") throw new Exception("Broker must run as SYSTEM");
            PipeSecurity security=new PipeSecurity(); security.SetAccessRuleProtection(true,false);
            security.SetOwner(new SecurityIdentifier("S-1-5-18"));
            security.AddAccessRule(new PipeAccessRule(new SecurityIdentifier("S-1-5-18"),PipeAccessRights.FullControl,AccessControlType.Allow));
            security.AddAccessRule(new PipeAccessRule(new SecurityIdentifier("S-1-5-32-544"),PipeAccessRights.FullControl,AccessControlType.Allow));
            security.AddAccessRule(new PipeAccessRule(new SecurityIdentifier(UserSid),PipeAccessRights.ReadWrite,AccessControlType.Allow));
            Audit(Reply("BrokerStarted"));
            while(true)
            {
                using(NamedPipeServerStream pipe=new NamedPipeServerStream(PipeName,PipeDirection.InOut,1,
                    PipeTransmissionMode.Byte,PipeOptions.Asynchronous,16384,16384,security))
                {
                    pipe.WaitForConnection();
                    try
                    {
                        int length=BitConverter.ToInt32(ReadBytes(pipe,4),0);
                        if(length<=0 || length>8192) throw new IOException("Invalid request length");
                        Dictionary<string,object> request=Json.Deserialize<Dictionary<string,object>>(Encoding.UTF8.GetString(ReadBytes(pipe,length)));
                        object result;
                        try { result=Handle(request); }
                        catch(Exception error) { Dictionary<string,object> failed=Reply("Failed"); failed["Error"]=error.Message; result=failed; }
                        Audit(new Dictionary<string,object> { {"Action",Text(request,"Action")}, {"Result",result} });
                        byte[] response=Encoding.UTF8.GetBytes(Json.Serialize(result));
                        WriteBytes(pipe,BitConverter.GetBytes(response.Length)); WriteBytes(pipe,response); pipe.Flush();
                    }
                    catch(Exception error) { Dictionary<string,object> rejected=Reply("RejectedConnection"); rejected["Error"]=error.Message; Audit(rejected); }
                }
            }
        }
        catch(Exception error) { try { Audit(new Dictionary<string,object>{{"Status","Fatal"},{"Error",error.ToString()}}); } catch {} return 1; }
    }
}
