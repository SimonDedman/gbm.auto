
UINT32 org.freedesktop.Notifications.Notify (STRING app_name, UINT32 replaces_id, STRING app_icon, STRING summary, STRING body, ARRAY actions, DICT hints, INT32 expire_timeout);
app_name	STRING	The optional name of the application sending the notification. Can be blank. 
summary	STRING	The summary text briefly describing the notification.
body	STRING	The optional detailed body text. Can be empty.
expire_timeout	INT32	The timeout time in milliseconds since the display of the notification at which the notification should automatically close.
If -1, the notification's expiration time is dependent on the notification server's settings, and may vary for the type of notification. If 0, never expire. 

UINT32 org.freedesktop.Notifications.Notify ("RStudio" app_name, "Processing Complete" summary, "Your code has finished processing" body, -1 expire_timeout);
system() in R
