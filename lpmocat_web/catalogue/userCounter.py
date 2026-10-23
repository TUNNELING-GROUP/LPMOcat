import time
from django.core.cache import cache

ACTIVE_USERS_KEY = "active_users"
ACTIVE_WINDOW = 10

class ActiveUsersMid:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        if not request.session.session_key:
            request.session.create()

        session_key = request.session.session_key
        now = time.time()

        active_users = cache.get("active_users", {})
        active_users[session_key] = now

        active_users = {
            k: v for k, v in active_users.items()
            if now - v < ACTIVE_WINDOW
        }

        cache.set("active_users", active_users, timeout=ACTIVE_WINDOW)

        request.active_users_count = len(active_users)
        # print("ACTIVE SESSIONS: ", request.session.session_key)

        return self.get_response(request)
